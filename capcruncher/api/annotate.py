import warnings
from typing import Literal, Union

import numpy as np
import pandas as pd
import pyranges1 as pr
from pandas.api.types import is_numeric_dtype

from capcruncher.utils import convert_bed_to_dataframe, convert_bed_to_pr

warnings.simplefilter("ignore", category=RuntimeWarning)


ROW_ID_COLUMN = "__cc_row_id"
INTERVAL_COLUMNS = ["Chromosome", "Start", "End", "Name"]
ANNOTATION_EXCLUDE_COLUMNS = {
    *INTERVAL_COLUMNS,
    ROW_ID_COLUMN,
}


def _as_pyranges(bed: Union[str, pd.DataFrame, pr.PyRanges]) -> pr.PyRanges:
    """Normalize supported BED-like inputs to a PyRanges1 interval frame.

    PyRanges1 objects are pandas DataFrame subclasses, so downstream pandas
    operations should work on the PyRanges frame directly rather than coercing
    it back out through ``pd.DataFrame(...)``.
    """

    if isinstance(bed, pr.PyRanges):
        return bed.copy()

    return convert_bed_to_pr(bed)


def _annotation_name_dtype(annotations: pr.PyRanges) -> object:
    """Return a nullable dtype suitable for values copied from annotation Name."""

    if "Name" not in annotations.columns:
        return pd.StringDtype()

    names = annotations["Name"]
    nullable_numeric_dtypes = {
        "int64": "Int64",
        "int32": "Int32",
        "float64": "Float64",
        "float32": "Float32",
    }
    if is_numeric_dtype(names):
        return nullable_numeric_dtypes.get(str(names.dtype), names.dtype)

    if isinstance(names.dtype, pd.CategoricalDtype):
        if is_numeric_dtype(names.cat.categories):
            return nullable_numeric_dtypes.get(str(names.cat.categories.dtype), "Int64")
        return names.dtype

    return pd.CategoricalDtype([*names.dropna().unique().astype(str)])


def _add_row_ids(intervals: pr.PyRanges) -> pr.PyRanges:
    return intervals.copy().assign(**{ROW_ID_COLUMN: np.arange(intervals.shape[0])})


def _prepare_annotation_intervals(
    annotations: Union[str, pd.DataFrame, pr.PyRanges],
    fallback_name: str,
) -> pr.PyRanges:
    intervals = _as_pyranges(annotations)
    if intervals.empty:
        return intervals

    intervals = intervals.copy()
    if "Name" not in intervals.columns:
        intervals["Name"] = [
            f"{fallback_name}_{idx}" for idx in range(intervals.shape[0])
        ]

    return intervals


def _split_query_metadata(
    query: pr.PyRanges,
) -> tuple[pr.PyRanges, dict[int, object], pd.DataFrame]:
    query = _add_row_ids(query)
    original_names = dict(zip(query[ROW_ID_COLUMN], query["Name"]))
    metadata_columns = [
        column for column in query.columns if column not in ANNOTATION_EXCLUDE_COLUMNS
    ]
    metadata = (
        query.set_index(ROW_ID_COLUMN).loc[:, metadata_columns]
        if metadata_columns
        else pd.DataFrame()
    )

    return (
        query.loc[:, [*INTERVAL_COLUMNS, ROW_ID_COLUMN]],
        original_names,
        metadata,
    )


def _restore_query_metadata(
    annotated: pr.PyRanges,
    original_names: dict[int, object],
    metadata: pd.DataFrame,
) -> pr.PyRanges:
    if not metadata.empty:
        annotated = (
            annotated.set_index(ROW_ID_COLUMN)
            .join(metadata, how="left")
            .reset_index()
        )

    return (
        annotated.assign(Name=lambda df: df[ROW_ID_COLUMN].map(original_names))
        .drop(columns=ROW_ID_COLUMN, errors="ignore")
        .reset_index(drop=True)
    )


def _overlaps(
    query: pr.PyRanges,
    annotations: pr.PyRanges,
    fraction: float,
) -> pr.PyRanges:
    if annotations.empty:
        return annotations

    overlaps = query.join_overlaps(
        annotations,
        strand_behavior="ignore",
        report_overlap_column="Overlap",
        suffix="_b",
    )
    if overlaps.empty:
        return overlaps

    return overlaps.assign(
        FractionOverlap=lambda df: df["Overlap"] / (df["End"] - df["Start"])
    ).loc[lambda df: df["FractionOverlap"] >= fraction]


def _get_annotation(
    query: pr.PyRanges,
    annotations: pr.PyRanges,
    name: str,
    fraction: float,
) -> pr.PyRanges:
    dtype = _annotation_name_dtype(annotations)
    hits = _overlaps(query, annotations, fraction)

    if hits.empty:
        values = pd.Series(dtype=pd.StringDtype(), name=name)
    else:
        values = (
            hits.drop_duplicates(subset=ROW_ID_COLUMN, keep="first")
            .set_index(ROW_ID_COLUMN)["Name_b"]
            .rename(name)
        )

    result = query.copy()
    result[name] = result[ROW_ID_COLUMN].map(values)
    result[name] = result[name].astype(dtype)
    return result


def _count_annotations(
    query: pr.PyRanges,
    annotations: pr.PyRanges,
    name: str,
    fraction: float,
) -> pr.PyRanges:
    hits = _overlaps(query, annotations, fraction)
    if hits.empty:
        counts = pd.Series(dtype=pd.Int8Dtype(), name=name)
    else:
        counts = hits.groupby(ROW_ID_COLUMN).size().astype(pd.Int8Dtype()).rename(name)

    return query.assign(
        **{name: lambda df: df[ROW_ID_COLUMN].map(counts).fillna(0).astype("Int8")}
    )


def _failed_annotation(query: pr.PyRanges, name: str) -> pr.PyRanges:
    return (
        query.assign(**{name: pd.NA})
        .assign(**{name: lambda df: df[name].astype(pd.StringDtype())})
    )


def annotate_intervals(
    query: Union[str, pd.DataFrame, pr.PyRanges],
    annotations: Union[str, pd.DataFrame, pr.PyRanges],
    name: str,
    method: Literal["get", "count"] = "get",
    fraction: float = 0,
    tolerate_errors: bool = True,
) -> pr.PyRanges:
    """Annotate a BED-like interval table with another BED-like table.

    The returned frame has one row per input query interval, preserves query
    row order and metadata columns, and stores either annotation names (`get`)
    or annotation counts (`count`) in ``name``.
    """

    prepared_query, original_names, metadata = _split_query_metadata(
        _as_pyranges(query)
    )

    try:
        annotation_intervals = _prepare_annotation_intervals(
            annotations, fallback_name=name
        )
        if method == "get":
            annotated = _get_annotation(
                prepared_query, annotation_intervals, name, fraction
            )
        elif method == "count":
            annotated = _count_annotations(
                prepared_query, annotation_intervals, name, fraction
            )
        else:
            annotated = _failed_annotation(prepared_query, name)
    except (
        OSError,
        IndexError,
        FileNotFoundError,
        StopIteration,
        AssertionError,
        TypeError,
        ValueError,
    ):
        if not tolerate_errors:
            raise
        annotated = _failed_annotation(prepared_query, name)

    return _restore_query_metadata(annotated, original_names, metadata)


def increase_cis_slice_priority(df: pd.DataFrame, score_multiplier: float = 2):
    """Prioritize cis slices by increasing their mapping score."""

    df = df.copy()
    df["parent_name"] = df["name"].str.split("|").str[0]

    df_chrom_counts = (
        df[["parent_name", "chrom"]].value_counts().to_frame("slices_per_chrom")
    )
    modal_chrom = (
        df_chrom_counts.groupby("parent_name")["slices_per_chrom"]
        .transform("max")
        .reset_index()
        .set_index("parent_name")["chrom"]
        .to_dict()
    )
    df["fragment_chrom"] = df["parent_name"].map(modal_chrom)
    df["score"] = np.where(
        df["chrom"] == df["fragment_chrom"],
        df["score"] * score_multiplier,
        df["score"] / score_multiplier,
    )

    return df.drop(columns="parent_name")


def remove_duplicates_from_bed(
    bed: pr.PyRanges,
    prioritize_cis_slices: bool = False,
    chroms_to_prioritize: Union[list, np.ndarray] = None,
) -> pr.PyRanges:
    """Remove duplicate BED names, using deterministic random tie-breaking."""

    df = convert_bed_to_dataframe(bed)
    if df.empty:
        return pr.PyRanges()

    if "name" not in df.columns:
        df = df.copy()
        df["name"] = [f"slice_{idx}" for idx in range(df.shape[0])]

    df = df.sample(frac=1, random_state=0)

    if prioritize_cis_slices:
        df = increase_cis_slice_priority(df)

    sort_columns = []
    sort_ascending = []
    if "score" in df.columns:
        sort_columns.append("score")
        sort_ascending.append(False)

    if chroms_to_prioritize:
        df = df.assign(
            is_chrom_priority=lambda frame: (
                frame["chrom"].isin(chroms_to_prioritize).astype(int)
            )
        )
        sort_columns.append("is_chrom_priority")
        sort_ascending.append(False)

    if sort_columns:
        df = df.sort_values(sort_columns, ascending=sort_ascending, kind="mergesort")

    df = (
        df.drop_duplicates(subset="name", keep="first")
        .sort_values(["chrom", "start"], kind="mergesort")
        .loc[:, ["chrom", "start", "end", "name"]]
        .rename(
            columns={
                "chrom": "Chromosome",
                "start": "Start",
                "end": "End",
                "name": "Name",
            }
        )
    )
    return pr.PyRanges(df)
