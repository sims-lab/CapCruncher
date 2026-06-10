"""CapCruncher utility library.

Submodules:
    bed       — BED file I/O, validation, and interval intersection
    genomics  — restriction enzymes, GTF and coordinate conversion
    io        — dict serialisation, cooler URIs, file-type detection
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Iterable
from functools import wraps
from typing import Any

import pandas as pd

from capcruncher.types import FLAG_NONE_VALUES, FLAG_OFF_VALUES, FLAG_ON_VALUES
from capcruncher.utils.bed import (
    BED_COLUMN_ALIASES,
    BED_COLUMN_CASE,
    BED_COLUMN_NAMES,
    INTERSECT_COLUMNS,
    BedInput,
    BedSchema,
    _prepare_intersection_frame,
    _read_bed_dataframe,
    _standardize_bed_columns,
    bed_has_duplicate_names,
    bed_has_name,
    convert_bed_to_dataframe,
    convert_bed_to_pr,
    format_coordinates,
    intersect_bins,
    is_valid_bed,
    split_intervals_on_chrom,
    validate_bed_dataframe,
)
from capcruncher.utils.genomics import (
    convert_interval_to_coords,
    get_human_readable_number_of_bp,
    get_restriction_site,
    gtf_line_to_bed12_line,
)
from capcruncher.utils.io import (
    get_cooler_uri,
    get_file_type,
    is_tabix,
    load_dict,
    read_dataframes,
    save_dict,
)

__all__ = [
    # bed
    "BED_COLUMN_ALIASES",
    "BED_COLUMN_CASE",
    "BED_COLUMN_NAMES",
    "INTERSECT_COLUMNS",
    "BedInput",
    "BedSchema",
    "_prepare_intersection_frame",
    "_read_bed_dataframe",
    "_standardize_bed_columns",
    "bed_has_duplicate_names",
    "bed_has_name",
    "convert_bed_to_dataframe",
    "convert_bed_to_pr",
    "format_coordinates",
    "intersect_bins",
    "is_valid_bed",
    "split_intervals_on_chrom",
    "validate_bed_dataframe",
    # genomics
    "convert_interval_to_coords",
    "get_human_readable_number_of_bp",
    "get_restriction_site",
    "gtf_line_to_bed12_line",
    # io
    "get_cooler_uri",
    "get_file_type",
    "is_tabix",
    "load_dict",
    "read_dataframes",
    "save_dict",
    # misc (defined below)
    "cycle_argument",
    "is_on",
    "is_off",
    "is_none",
    "hash_column",
    "get_timing",
    "categorise_tracks",
]


# ---------------------------------------------------------------------------
# Miscellaneous helpers
# ---------------------------------------------------------------------------


def cycle_argument(arg: list) -> Iterable:
    """Allows for the same argument to be stated once but repeated for all files"""

    if len(arg) == 1:
        return itertools.cycle((arg[0],))
    else:
        return arg


def is_on(param: str) -> bool:
    return str(param).strip().lower() in FLAG_ON_VALUES


def is_off(param: str) -> bool:
    return str(param).strip().lower() in FLAG_OFF_VALUES


def is_none(param: str) -> bool:
    return str(param).strip().lower() in FLAG_NONE_VALUES


def hash_column(col: Iterable, hash_type: int = 64) -> list:
    """Hashing using xxhash on an iterable. Not vectorised."""
    import xxhash

    hash_dict = {
        32: xxhash.xxh32_intdigest,
        64: xxhash.xxh64_intdigest,
        128: xxhash.xxh128_intdigest,
    }

    hash_func = hash_dict.get(hash_type)
    if hash_func is None:
        raise ValueError(f"Unsupported hash type: {hash_type}")

    return [hash_func(v) for v in col]


def get_timing(task_name: str | None = None) -> Callable:
    """Decorator: records the time taken by the wrapped function."""
    import time
    from datetime import timedelta

    from loguru import logger

    def wrapper(f: Callable) -> Callable:
        @wraps(f)
        def wrapped(*args: Any, **kwargs: Any) -> Any:
            time_start = time.perf_counter()
            result = f(*args, **kwargs)
            time_end = time.perf_counter()

            time_taken = timedelta(seconds=(time_end - time_start))
            logger.info(f"Completed {task_name} in {time_taken} (hh:mm:ss.ms)")
            return result

        return wrapped

    return wrapper


def categorise_tracks(ser: pd.Series) -> list:
    """Gets a series for grouping tracks together"""
    mapping = {
        "raw": "Replicates",
        "normalised": "Replicates_Scaled",
        "norm": "Replicates_Scaled",
        "summary": "Samples_Summarised",
        "subtraction": "Samples_Compared",
    }
    categories = []
    for _, value in ser.items():
        for key in mapping:
            if key in value:
                categories.append(mapping[key])

    return categories
