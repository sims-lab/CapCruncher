from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import pandas as pd
import polars as pl

from capcruncher.types import Assay


type FilterFunction = Callable[[pl.DataFrame], pl.DataFrame]


class FilterStepName(StrEnum):
    GET_UNFILTERED_SLICES = "get_unfiltered_slices"
    REMOVE_UNMAPPED_SLICES = "remove_unmapped_slices"
    REMOVE_ORPHAN_SLICES = "remove_orphan_slices"
    REMOVE_DUPLICATE_RE_FRAGS = "remove_duplicate_re_frags"
    REMOVE_SLICES_WITHOUT_RE_FRAG_ASSIGNED = "remove_slices_without_re_frag_assigned"
    REMOVE_DUPLICATE_SLICES = "remove_duplicate_slices"
    REMOVE_DUPLICATE_SLICES_PE = "remove_duplicate_slices_pe"
    REMOVE_EXCLUDED_SLICES = "remove_excluded_slices"
    REMOVE_BLACKLISTED_SLICES = "remove_blacklisted_slices"
    REMOVE_NON_REPORTER_FRAGMENTS = "remove_non_reporter_fragments"
    REMOVE_MULTI_CAPTURE_FRAGMENTS = "remove_multi_capture_fragments"
    REMOVE_VIEWPOINT_ADJACENT_RESTRICTION_FRAGMENTS = (
        "remove_viewpoint_adjacent_restriction_fragments"
    )
    REMOVE_SLICES_WITH_ONE_REPORTER = "remove_slices_with_one_reporter"
    REMOVE_SLICES_OUTSIDE_CAPTURE = "remove_slices_outside_capture"
    REMOVE_NON_CAPTURE_FRAGMENTS = "remove_non_capture_fragments"
    REMOVE_DUAL_CAPTURE_FRAGMENTS = "remove_dual_capture_fragments"
    REMOVE_RELIGATION = "remove_religation"

    @classmethod
    def parse(cls, value: str) -> "FilterStepName":
        try:
            return cls(value)
        except ValueError as exc:
            valid = ", ".join(step.value for step in cls)
            raise ValueError(
                f"Unknown filter step {value!r}. Valid steps are: {valid}."
            ) from exc


type FilterStepInput = str | FilterStepName


def _coerce_filter_step_name(value: str | FilterStepName) -> FilterStepName:
    if isinstance(value, FilterStepName):
        return value
    return FilterStepName.parse(value)


@dataclass(frozen=True)
class RegisteredFilterStep:
    name: FilterStepName
    function: FilterFunction
    assays: frozenset[Assay]


class FilterStepRegistry:
    def __init__(self) -> None:
        self._steps: dict[FilterStepName, RegisteredFilterStep] = {}

    def register(
        self,
        name: FilterStepName,
        function: FilterFunction,
        *,
        assays: Sequence[Assay] = tuple(Assay),
    ) -> None:
        self._steps[name] = RegisteredFilterStep(name, function, frozenset(assays))

    def validate_plan(self, plan: Any) -> None:
        for stage in plan.stages:
            for step_name in stage.steps:
                step = self._steps.get(step_name)
                if step is None:
                    valid = ", ".join(sorted(step.value for step in self._steps))
                    raise ValueError(
                        f"Unknown filter step {step_name.value!r}. "
                        f"Valid steps are: {valid}."
                    )
                if plan.assay not in step.assays:
                    valid_assays = ", ".join(sorted(assay.value for assay in step.assays))
                    raise ValueError(
                        f"Filter step {step_name.value!r} is not valid for assay "
                        f"{plan.assay.value!r}. Valid assays: {valid_assays}."
                    )

    def get(self, name: FilterStepName) -> FilterFunction:
        return self._steps[name].function

    @property
    def names(self) -> list[str]:
        return sorted(step.value for step in self._steps)


def _normalise_slices(slices: pd.DataFrame | pl.DataFrame) -> pl.DataFrame:
    df = slices.clone() if isinstance(slices, pl.DataFrame) else pl.from_pandas(slices)

    fill_zero_columns = [
        column
        for column in ("blacklist", "capture_count", "exclusion_count")
        if column in df.columns
    ]
    if fill_zero_columns:
        df = df.with_columns(pl.col(column).fill_null(0) for column in fill_zero_columns)

    casts: list[pl.Expr] = []
    if "blacklist" in df.columns:
        casts.append(pl.col("blacklist").cast(pl.Float64, strict=False))
    if "capture_count" in df.columns:
        casts.append(pl.col("capture_count").cast(pl.Float64, strict=False))
    if "exclusion_count" in df.columns:
        casts.append(pl.col("exclusion_count").cast(pl.Float64, strict=False))
    if "restriction_fragment" in df.columns:
        casts.append(pl.col("restriction_fragment").cast(pl.Int64, strict=False))
    if "mapped" in df.columns:
        casts.append(pl.col("mapped").cast(pl.Boolean, strict=False))
    if "multimapped" in df.columns:
        casts.append(pl.col("multimapped").cast(pl.Boolean, strict=False))
    if casts:
        df = df.with_columns(casts)

    sort_columns = [column for column in ("parent_read", "slice") if column in df.columns]
    return df.sort(sort_columns) if sort_columns else df


def _semi_join_parent_ids(df: pl.DataFrame, parent_ids: pl.DataFrame) -> pl.DataFrame:
    return df.join(parent_ids.select("parent_id").unique(), on="parent_id", how="semi")


def _anti_join_parent_ids(df: pl.DataFrame, parent_ids: pl.DataFrame) -> pl.DataFrame:
    return df.join(parent_ids.select("parent_id").unique(), on="parent_id", how="anti")


def _fragment_coordinate_signatures(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.sort(["parent_id", "slice"])
        .group_by("parent_id", maintain_order=True)
        .agg(pl.col("coordinates").cast(pl.Utf8).str.join("|").alias("coords"))
    )


def _captures(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter(pl.col("capture_count") == 1)


def _slices_with_viewpoint(df: pl.DataFrame) -> pl.DataFrame:
    captures = _captures(df).select(
        pl.col("parent_id"),
        pl.col("capture").alias("viewpoint"),
    )
    return df.join(captures, on="parent_id", how="inner")


def get_unfiltered_slices(df: pl.DataFrame) -> pl.DataFrame:
    return df


def remove_unmapped_slices(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter(pl.col("mapped"))


def remove_orphan_slices(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.with_columns(pl.len().over("parent_id").alias("__parent_slice_count"))
        .filter(pl.col("__parent_slice_count") > 1)
        .drop("__parent_slice_count")
    )


def remove_duplicate_re_frags(df: pl.DataFrame) -> pl.DataFrame:
    return df.unique(
        subset=["parent_read", "restriction_fragment"],
        keep="first",
        maintain_order=True,
    )


def remove_slices_without_re_frag_assigned(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter(pl.col("restriction_fragment").is_not_null())


def remove_duplicate_slices(df: pl.DataFrame) -> pl.DataFrame:
    deduplicated = _fragment_coordinate_signatures(df).unique(
        subset=["coords"], keep="first", maintain_order=True
    )
    return _semi_join_parent_ids(df, deduplicated)


def remove_duplicate_slices_pe(df: pl.DataFrame) -> pl.DataFrame:
    if df.is_empty() or "pe" not in df.columns:
        return df

    has_pe = (
        df.head(100)
        .select(pl.col("pe").cast(pl.Utf8).str.contains("pe").sum())
        .item()
    )
    if has_pe <= 1:
        return df

    fragments = _fragment_coordinate_signatures(df).with_columns(
        read_start=pl.col("coords")
        .str.split("|")
        .list.first()
        .str.extract(r":(\d+)-", 1),
        read_end=pl.col("coords")
        .str.split("|")
        .list.last()
        .str.extract(r"-(\d+)$", 1),
    )
    deduplicated = fragments.unique(
        subset=["read_start", "read_end"], keep="first", maintain_order=True
    )
    return _semi_join_parent_ids(df, deduplicated)


def remove_excluded_slices(df: pl.DataFrame) -> pl.DataFrame:
    with_viewpoint = _slices_with_viewpoint(df)
    passed = with_viewpoint.filter(
        (pl.col("exclusion_count") < 1)
        | (
            pl.col("exclusion").cast(pl.Utf8).fill_null("")
            != pl.col("viewpoint").cast(pl.Utf8).fill_null("")
        )
    )
    return _semi_join_parent_ids(df, passed)


def remove_blacklisted_slices(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter((pl.col("blacklist") == 0) | pl.col("blacklist").is_null())


def remove_non_reporter_fragments(df: pl.DataFrame) -> pl.DataFrame:
    fragments = df.group_by("parent_id", maintain_order=True).agg(
        n_capture=pl.col("capture_count").sum(),
        n_mapped=pl.col("mapped").cast(pl.Int64).sum(),
        n_blacklist=pl.col("blacklist").sum(),
        n_exclusions=pl.col("exclusion_count").sum(),
    )
    with_reporters = fragments.filter(
        (pl.col("n_mapped") - pl.col("n_capture") - pl.col("n_blacklist") - pl.col("n_exclusions"))
        > 0
    )
    return _semi_join_parent_ids(df, with_reporters)


def remove_multi_capture_fragments(df: pl.DataFrame) -> pl.DataFrame:
    single_capture = (
        df.filter(pl.col("capture").is_not_null())
        .group_by("parent_id", maintain_order=True)
        .agg(pl.col("capture").n_unique().alias("n_captures"))
        .filter(pl.col("n_captures") == 1)
    )
    return _semi_join_parent_ids(df, single_capture)


def remove_viewpoint_adjacent_restriction_fragments(
    df: pl.DataFrame, n_adjacent: int = 1
) -> pl.DataFrame:
    captures = (
        _captures(df)
        .select("capture", "restriction_fragment")
        .unique(subset=["capture"], keep="first", maintain_order=True)
        .with_columns(
            exclusion_start=pl.col("restriction_fragment") - n_adjacent,
            exclusion_end=pl.col("restriction_fragment") + n_adjacent,
        )
    )
    if captures.is_empty():
        return df

    excluded = (
        _slices_with_viewpoint(df)
        .join(captures, left_on="viewpoint", right_on="capture", how="inner")
        .filter(
            (pl.col("restriction_fragment") >= pl.col("exclusion_start"))
            & (pl.col("restriction_fragment") <= pl.col("exclusion_end"))
            & (pl.col("capture_count") == 0)
        )
    )
    return _anti_join_parent_ids(df, excluded)


def _capture_fragments(df: pl.DataFrame) -> pl.DataFrame:
    fragments = (
        df.sort(["parent_read", "chrom", "start"])
        .group_by("parent_read", maintain_order=True)
        .agg(
            pl.col("slice").n_unique().alias("unique_slices"),
            pl.col("pe").first().alias("pe"),
            pl.col("mapped").cast(pl.Int64).sum().alias("mapped"),
            pl.col("multimapped").cast(pl.Int64).sum().alias("multimapped"),
            pl.col("capture").n_unique().alias("unique_capture_sites"),
            pl.col("capture_count").sum().alias("capture_count"),
            pl.col("exclusion").n_unique().alias("unique_exclusions"),
            pl.col("exclusion_count").sum().alias("exclusion_count"),
            pl.col("restriction_fragment").n_unique().alias("unique_restriction_fragments"),
            pl.col("blacklist").sum().alias("blacklist"),
            pl.col("coordinates").cast(pl.Utf8).str.join("|").alias("coordinates"),
        )
        .with_columns(
            reporter_count=pl.when(pl.col("capture_count") > 0)
            .then(
                pl.col("mapped")
                - (
                    pl.col("exclusion_count")
                    + pl.col("capture_count")
                    + pl.col("blacklist")
                )
            )
            .otherwise(0)
        )
    )
    return fragments


def _tiled_fragments(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.sort(["parent_read", "chrom", "start"])
        .group_by("parent_read", maintain_order=True)
        .agg(
            pl.col("parent_id").first().alias("id"),
            pl.col("slice").n_unique().alias("unique_slices"),
            pl.col("pe").first().alias("pe"),
            pl.col("mapped").cast(pl.Int64).sum().alias("mapped"),
            pl.col("multimapped").cast(pl.Int64).sum().alias("multimapped"),
            pl.col("capture_count").sum().alias("capture_count"),
            pl.col("restriction_fragment").n_unique().alias("unique_restriction_fragments"),
            pl.col("blacklist").sum().alias("blacklisted_slices"),
            pl.col("coordinates").cast(pl.Utf8).str.join("|").alias("coordinates"),
        )
    )


def remove_slices_with_one_reporter(df: pl.DataFrame) -> pl.DataFrame:
    fragments = _capture_fragments(df).filter(pl.col("reporter_count") > 1)
    return df.join(
        fragments.select("parent_read").unique(), on="parent_read", how="semi"
    )


def remove_slices_outside_capture(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter(pl.col("capture_count") != 0)


def remove_non_capture_fragments(df: pl.DataFrame) -> pl.DataFrame:
    with_capture = (
        df.group_by("parent_id", maintain_order=True)
        .agg(pl.col("capture_count").sum().alias("capture_count"))
        .filter(pl.col("capture_count") > 0)
    )
    return _semi_join_parent_ids(df, with_capture)


def remove_dual_capture_fragments(df: pl.DataFrame) -> pl.DataFrame:
    single_capture = (
        _captures(df)
        .group_by("parent_id", maintain_order=True)
        .agg(pl.col("capture").n_unique().alias("n_captures"))
        .filter(pl.col("n_captures") <= 1)
    )
    return _semi_join_parent_ids(df, single_capture)


def remove_religation(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.sort("restriction_fragment")
        .with_columns(
            pl.col("restriction_fragment")
            .diff()
            .over("parent_id")
            .fill_null(-1)
            .alias("__restriction_fragment_diff")
        )
        .filter(
            (pl.col("__restriction_fragment_diff") > 1)
            | (pl.col("__restriction_fragment_diff") == -1)
        )
        .drop("__restriction_fragment_diff")
    )


FILTER_REGISTRY = FilterStepRegistry()
for _name, _function in {
    FilterStepName.GET_UNFILTERED_SLICES: get_unfiltered_slices,
    FilterStepName.REMOVE_UNMAPPED_SLICES: remove_unmapped_slices,
    FilterStepName.REMOVE_ORPHAN_SLICES: remove_orphan_slices,
    FilterStepName.REMOVE_DUPLICATE_RE_FRAGS: remove_duplicate_re_frags,
    FilterStepName.REMOVE_SLICES_WITHOUT_RE_FRAG_ASSIGNED: (
        remove_slices_without_re_frag_assigned
    ),
    FilterStepName.REMOVE_DUPLICATE_SLICES: remove_duplicate_slices,
    FilterStepName.REMOVE_DUPLICATE_SLICES_PE: remove_duplicate_slices_pe,
    FilterStepName.REMOVE_BLACKLISTED_SLICES: remove_blacklisted_slices,
}.items():
    FILTER_REGISTRY.register(_name, _function)

FILTER_REGISTRY.register(
    FilterStepName.REMOVE_EXCLUDED_SLICES,
    remove_excluded_slices,
    assays=(Assay.CAPTURE,),
)
FILTER_REGISTRY.register(
    FilterStepName.REMOVE_NON_REPORTER_FRAGMENTS,
    remove_non_reporter_fragments,
    assays=(Assay.CAPTURE, Assay.TRI),
)
FILTER_REGISTRY.register(
    FilterStepName.REMOVE_MULTI_CAPTURE_FRAGMENTS,
    remove_multi_capture_fragments,
    assays=(Assay.CAPTURE, Assay.TRI),
)
FILTER_REGISTRY.register(
    FilterStepName.REMOVE_VIEWPOINT_ADJACENT_RESTRICTION_FRAGMENTS,
    remove_viewpoint_adjacent_restriction_fragments,
    assays=(Assay.CAPTURE,),
)
FILTER_REGISTRY.register(
    FilterStepName.REMOVE_SLICES_WITH_ONE_REPORTER,
    remove_slices_with_one_reporter,
    assays=(Assay.TRI,),
)
FILTER_REGISTRY.register(
    FilterStepName.REMOVE_SLICES_OUTSIDE_CAPTURE,
    remove_slices_outside_capture,
    assays=(Assay.TILED,),
)
FILTER_REGISTRY.register(
    FilterStepName.REMOVE_NON_CAPTURE_FRAGMENTS,
    remove_non_capture_fragments,
    assays=(Assay.TILED,),
)
FILTER_REGISTRY.register(
    FilterStepName.REMOVE_DUAL_CAPTURE_FRAGMENTS,
    remove_dual_capture_fragments,
    assays=(Assay.TILED,),
)
FILTER_REGISTRY.register(
    FilterStepName.REMOVE_RELIGATION,
    remove_religation,
    assays=(Assay.TILED,),
)
