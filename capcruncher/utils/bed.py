"""BED file reading, validation, coordinate manipulation, and interval intersection."""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd
import pandera.pandas as pa
import pyranges1 as pr
from pandera import errors as pandera_errors
from pandera.typing.pandas import Series as PASeries


class BedSchema(pa.DataFrameModel):
    """Pandera schema for a minimal BED DataFrame (BED3+).

    Enforces chrom as str (avoids mixed int/str dtype — root cause of issue #234),
    non-negative coordinates, and end > start.
    """

    chrom: PASeries[str] = pa.Field(nullable=False)
    start: PASeries[int] = pa.Field(ge=0)
    end: PASeries[int] = pa.Field(ge=1)

    class Config:
        coerce = True
        strict = False  # extra columns (name, score, strand …) are allowed

    @pa.dataframe_check
    @classmethod
    def end_gt_start(cls, df: pd.DataFrame) -> pd.Series:
        return df["end"] > df["start"]


def validate_bed_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Validate and coerce a BED DataFrame against BedSchema.

    Returns the coerced DataFrame. Raises SchemaError on invalid data.
    Only validates when the required columns (chrom, start, end) are present.
    """
    if not {"chrom", "start", "end"}.issubset(df.columns):
        return df
    return BedSchema.validate(df)


type BedInput = str | os.PathLike | pd.DataFrame | pr.PyRanges

BED_COLUMN_NAMES = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "item_rgb",
    "block_count",
    "block_sizes",
    "block_starts",
]

BED_COLUMN_CASE = {
    "chrom": "Chromosome",
    "start": "Start",
    "end": "End",
    "name": "Name",
    "score": "Score",
    "strand": "Strand",
    "thick_start": "ThickStart",
    "thick_end": "ThickEnd",
    "item_rgb": "ItemRGB",
    "block_count": "BlockCount",
    "block_sizes": "BlockSizes",
    "block_starts": "BlockStarts",
}

BED_COLUMN_ALIASES = {
    "chrom": "chrom",
    "chromosome": "chrom",
    "start": "start",
    "end": "end",
    "name": "name",
    "score": "score",
    "strand": "strand",
    "thickstart": "thick_start",
    "thickend": "thick_end",
    "itemrgb": "item_rgb",
    "blockcount": "block_count",
    "blocksizes": "block_sizes",
    "blockstarts": "block_starts",
}

INTERSECT_COLUMNS = [
    "chrom_1",
    "start_1",
    "end_1",
    "name_1",
    "chrom_2",
    "start_2",
    "end_2",
    "name_2",
    "overlap",
]


def _read_bed_dataframe(bed: BedInput, nrows: int | None = None) -> pd.DataFrame:
    if isinstance(bed, pr.PyRanges):
        return bed.copy()

    if isinstance(bed, pd.DataFrame):
        return bed.copy()

    df = pd.read_csv(
        bed,
        sep="\t",
        header=None,
        comment="#",
        nrows=nrows,
        dtype={0: str},
    )

    df.columns = [
        BED_COLUMN_NAMES[i] if i < len(BED_COLUMN_NAMES) else f"col_{i}"
        for i in range(df.shape[1])
    ]
    return df


def _standardize_bed_columns(
    df: pd.DataFrame, capitalized: bool = False
) -> pd.DataFrame:
    rename_map = {}
    for column in df.columns:
        alias_key = re.sub(r"[^a-z0-9]", "", str(column).lower())
        canonical = BED_COLUMN_ALIASES.get(alias_key)
        if canonical:
            rename_map[column] = (
                BED_COLUMN_CASE[canonical] if capitalized else canonical
            )

    return df.rename(columns=rename_map)


def _prepare_intersection_frame(df: BedInput, name_prefix: str) -> pd.DataFrame:
    frame = convert_bed_to_dataframe(df)
    if frame.empty:
        return frame

    frame = _standardize_bed_columns(frame, capitalized=False)
    for column in ("start", "end"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="raise")
    if "name" not in frame.columns:
        frame = frame.copy()
        frame["name"] = [f"{name_prefix}_{idx}" for idx in range(frame.shape[0])]
    else:
        frame["name"] = frame["name"].fillna(
            pd.Series(
                [f"{name_prefix}_{idx}" for idx in range(frame.shape[0])],
                index=frame.index,
            )
        )

    return frame


def is_valid_bed(bed: BedInput, verbose: bool = True) -> bool:
    from loguru import logger

    """Return True when the first non-empty row has at least three BED columns."""

    try:
        df = _read_bed_dataframe(bed, nrows=1)
    except FileNotFoundError:
        if verbose:
            logger.warning(f"Bed file: {bed} not found")
        return False
    except pd.errors.EmptyDataError:
        if verbose:
            logger.warning(f"Bed file: {bed} is empty")
        return False
    except Exception as e:
        if verbose:
            logger.warning(f"Exception raised {e}")
        return False

    return df.shape[1] >= 3


def bed_has_name(bed: BedInput) -> bool:
    """Return True when the first non-empty row has at least four BED columns."""

    try:
        df = _read_bed_dataframe(bed, nrows=1)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return False

    return df.shape[1] >= 4


def bed_has_duplicate_names(bed: BedInput) -> bool:
    """Return True when a BED-like input has duplicate name values."""

    df = convert_bed_to_dataframe(bed)
    if "name" not in df.columns or df.empty:
        return False

    return df["name"].dropna().duplicated().any()


def split_intervals_on_chrom(intervals: BedInput) -> dict:
    """Creates dictionary from bed file with the chroms as keys"""

    intervals = convert_bed_to_dataframe(intervals)
    if intervals.empty or "chrom" not in intervals.columns:
        return {}

    return {chrom: df for chrom, df in intervals.groupby("chrom")}


def intersect_bins(
    bins_1: pd.DataFrame, bins_2: pd.DataFrame, **bedtools_kwargs
) -> pd.DataFrame:
    """Intersect two interval tables and return a labeled pandas DataFrame."""

    left = _prepare_intersection_frame(bins_1, name_prefix="region_1")
    right = _prepare_intersection_frame(bins_2, name_prefix="region_2")

    if left.empty or right.empty:
        return pd.DataFrame(columns=INTERSECT_COLUMNS)

    left = left.copy()
    right = right.copy()

    slack = int(bedtools_kwargs.get("slack", 0) or 0)
    if slack:
        left["start"] = (left["start"] - slack).clip(lower=0)
        left["end"] = left["end"] + slack
        right["start"] = (right["start"] - slack).clip(lower=0)
        right["end"] = right["end"] + slack

    strandedness = bedtools_kwargs.get("strandedness")
    if bedtools_kwargs.get("s"):
        strandedness = "same"

    joined = convert_bed_to_pr(left).join_overlaps(
        convert_bed_to_pr(right),
        strand_behavior="ignore",
        suffix="_2",
        report_overlap_column="overlap",
    )
    df_intersect = joined.copy()
    if df_intersect.empty:
        return pd.DataFrame(columns=INTERSECT_COLUMNS)

    if strandedness in {"same", "opposite"} and {"Strand", "Strand_2"}.issubset(
        df_intersect.columns
    ):
        if strandedness == "same":
            df_intersect = df_intersect[
                df_intersect["Strand"] == df_intersect["Strand_2"]
            ]
        else:
            df_intersect = df_intersect[
                df_intersect["Strand"] != df_intersect["Strand_2"]
            ]
        if df_intersect.empty:
            return pd.DataFrame(columns=INTERSECT_COLUMNS)

    return pd.DataFrame(
        {
            "chrom_1": df_intersect["Chromosome"],
            "start_1": df_intersect["Start"],
            "end_1": df_intersect["End"],
            "name_1": df_intersect["Name"],
            "chrom_2": df_intersect["Chromosome"],
            "start_2": df_intersect["Start_2"],
            "end_2": df_intersect["End_2"],
            "name_2": df_intersect["Name_2"],
            "overlap": df_intersect["overlap"],
        },
        columns=INTERSECT_COLUMNS,
    )


def convert_bed_to_pr(bed: BedInput) -> pr.PyRanges:
    """Convert a BED-like object to a PyRanges object."""

    df = convert_bed_to_dataframe(bed)
    if df.empty:
        return pr.PyRanges()

    df = _standardize_bed_columns(df, capitalized=True)
    for column in ("Start", "End"):
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="raise")
    if "Name" in df.columns:
        df["Name"] = df["Name"].astype("category")

    return pr.PyRanges(df)


def convert_bed_to_dataframe(bed: BedInput) -> pd.DataFrame:
    """Converts a BED-like object to a DataFrame-style interval table."""
    from loguru import logger

    if isinstance(bed, (str, os.PathLike)):
        try:
            bed_conv = _read_bed_dataframe(bed)
        except FileNotFoundError:
            logger.warning(f"File {bed} not found")
            bed_conv = pd.DataFrame()
        except pd.errors.EmptyDataError:
            logger.warning(f"File {bed} is empty")
            bed_conv = pd.DataFrame()

    elif isinstance(bed, pr.PyRanges):
        bed_conv = pd.DataFrame(bed)

    elif isinstance(bed, pd.DataFrame):
        bed_conv = bed.copy()

    else:
        raise TypeError(f"Unsupported BED input type: {type(bed)!r}")

    bed_conv = _standardize_bed_columns(bed_conv, capitalized=False)

    if not bed_conv.empty:
        try:
            bed_conv = validate_bed_dataframe(bed_conv)
        except pandera_errors.SchemaError:
            from loguru import logger

            if {"start", "end"}.issubset(bed_conv.columns):
                end_num = pd.Series(pd.to_numeric(bed_conv["end"], errors="coerce"))
                start_num = pd.Series(pd.to_numeric(bed_conv["start"], errors="coerce"))
                valid_mask = end_num.gt(start_num)
                n_dropped = int(len(valid_mask) - valid_mask.sum())
                if n_dropped:
                    logger.warning(f"Dropped {n_dropped} BED rows where end <= start")
                bed_conv = pd.DataFrame(bed_conv.loc[valid_mask]).reset_index(drop=True)
                if not bed_conv.empty:
                    bed_conv = validate_bed_dataframe(bed_conv)
            else:
                bed_conv = pd.DataFrame()

    return bed_conv


def format_coordinates(coordinates: str | os.PathLike) -> pr.PyRanges:
    """Convert coordinates supplied in string format or a BED file to PyRanges."""

    coordinates = str(coordinates)
    pattern_genomic_coord = re.compile(
        r"^(chr[0-2xXyYmM][0-9]*):(\d+)-(\d+)(?:\s+(\S+))?$"
    )

    match = pattern_genomic_coord.match(coordinates)
    if match:
        chrom, start, end, name = match.groups()
        if not name:
            name = "region_0"

        return pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": [chrom],
                    "Start": [int(start)],
                    "End": [int(end)],
                    "Name": [name],
                }
            )
        )

    path_name = Path(coordinates).name.lower()
    if path_name.endswith((".bed", ".bed.gz", ".bed.bgz")):
        if is_valid_bed(coordinates):
            bed_df = convert_bed_to_dataframe(coordinates)
            if bed_has_name(bed_df):
                return convert_bed_to_pr(bed_df)

            bed_df = bed_df[["chrom", "start", "end"]].copy()
            bed_df = bed_df.reset_index(drop=True)
            bed_df["name"] = bed_df.index.map(lambda idx: f"region_{idx}")
            return convert_bed_to_pr(bed_df)

        raise ValueError("Invalid bed file supplied.")

    raise ValueError(
        """Provide coordinates in the form chr[NUMBER]:[START]-[END]/BED file"""
    )
