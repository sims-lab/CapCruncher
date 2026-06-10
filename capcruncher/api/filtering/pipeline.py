from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import ClassVar

import pandas as pd
import polars as pl
from loguru import logger

from capcruncher.api.filtering.plan import DEFAULT_FILTER_PLANS, FilterPlan, FilterStage
from capcruncher.api.filtering.steps import (
    FILTER_REGISTRY,
    FilterStepInput,
    FilterStepName,
    FilterStepRegistry,
    _capture_fragments,
    _captures,
    _coerce_filter_step_name,
    _normalise_slices,
    _slices_with_viewpoint,
    _tiled_fragments,
    remove_viewpoint_adjacent_restriction_fragments,
)
from capcruncher.api.statistics import SliceFilterStats
from capcruncher.types import Assay, ReadType


class FilterPipeline:
    def __init__(
        self,
        slices: pd.DataFrame | pl.DataFrame,
        plan: FilterPlan,
        *,
        sample_name: str = "",
        read_type: ReadType | str = ReadType.FLASHED,
        registry: FilterStepRegistry = FILTER_REGISTRY,
    ) -> None:
        registry.validate_plan(plan)
        self.plan = plan
        self.registry = registry
        self.slices = _normalise_slices(slices)
        self.sample_name = sample_name
        self.read_type = (
            read_type.value if isinstance(read_type, ReadType) else str(read_type)
        )
        self.filtering_stats: list[SliceFilterStats] = []
        self.current_stage = ""
        self.current_filter: FilterStepName | None = None

    def run(
        self,
        *,
        output_slices: bool | str = False,
        output_location: Path | str = ".",
    ) -> None:
        output_path = Path(output_location)
        for stage in self.plan.stages:
            self.current_stage = stage.name
            for step_name in stage.steps:
                self.current_filter = step_name
                logger.info(f"Filtering slices: {step_name.value}")
                self.slices = self.registry.get(step_name)(self.slices)
                logger.info(f"Completed: {step_name.value}")
                logger.info(f"Number of slices: {self.slices.height}")
                logger.info(
                    "Number of reads: "
                    f"{self.slices.select(pl.col('parent_read').n_unique()).item()}"
                )

                if output_slices == "filter":
                    self.slices.write_csv(output_path / f"{step_name.value}.tsv")

            if output_slices == "stage":
                self.slices.write_csv(output_path / f"{stage.name}.tsv")

            self.filtering_stats.append(self.slice_stats(stage.name))

    def slice_stats(self, stage: str) -> SliceFilterStats:
        if self.slices.is_empty():
            n_slices = 0
            n_fragments = 0
        else:
            n_slices = self.slices.select(pl.col("slice_name").n_unique()).item()
            n_fragments = self.slices.select(pl.col("parent_read").n_unique()).item()

        return SliceFilterStats(
            sample=self.sample_name,
            stage=stage,
            n_fragments=n_fragments,
            n_slices=n_slices,
            read_type=self.read_type,
        )


class SliceFilter:
    assay: ClassVar[Assay]

    def __init__(
        self,
        slices: pd.DataFrame | pl.DataFrame,
        filter_stages: FilterPlan
        | Mapping[str, Sequence[FilterStepInput]]
        | None = None,
        sample_name: str = "",
        read_type: ReadType | str = ReadType.FLASHED,
        filter_profile: Path | str | None = None,
    ) -> None:
        self.sample_name = sample_name
        self.read_type = (
            read_type.value
            if isinstance(read_type, ReadType)
            else str(read_type or ReadType.FLASHED.value)
        )
        self.current_filter: FilterStepName | None = None
        self.current_stage = ""
        self.plan = self._resolve_filter_plan(filter_stages, filter_profile)
        self.pipeline = FilterPipeline(
            slices,
            self.plan,
            sample_name=sample_name,
            read_type=self.read_type,
        )
        self.filtering_stats: list[SliceFilterStats] = []

    def _resolve_filter_plan(
        self,
        filter_stages: FilterPlan | Mapping[str, Sequence[FilterStepInput]] | None,
        filter_profile: Path | str | None,
    ) -> FilterPlan:
        if filter_profile is not None:
            return FilterPlan.from_toml(filter_profile, expected_assay=self.assay)
        if isinstance(filter_stages, FilterPlan):
            if filter_stages.assay != self.assay:
                raise ValueError(
                    f"Filter plan assay {filter_stages.assay.value!r} does not match "
                    f"{self.assay.value!r}."
                )
            return filter_stages
        if isinstance(filter_stages, Mapping):
            return FilterPlan(
                assay=self.assay,
                stages=tuple(
                    FilterStage(
                        str(stage),
                        tuple(_coerce_filter_step_name(step) for step in steps),
                    )
                    for stage, steps in filter_stages.items()
                ),
            )
        if filter_stages is not None:
            return FilterPlan.from_toml(filter_stages, expected_assay=self.assay)
        return DEFAULT_FILTER_PLANS[self.assay]

    @property
    def _polars_slices(self) -> pl.DataFrame:
        return self.pipeline.slices

    @_polars_slices.setter
    def _polars_slices(self, slices: pl.DataFrame) -> None:
        self.pipeline.slices = slices

    @property
    def slices(self) -> pl.DataFrame:
        return self.pipeline.slices

    @slices.setter
    def slices(self, slices: pd.DataFrame | pl.DataFrame) -> None:
        self.pipeline.slices = _normalise_slices(slices)

    @property
    def filters(self) -> list[str]:
        return FILTER_REGISTRY.names

    @property
    def filter_stages(self) -> dict[str, list[str]]:
        return {
            stage.name: [step.value for step in stage.steps]
            for stage in self.plan.stages
        }

    @property
    def slice_stats(self) -> SliceFilterStats:
        return self.pipeline.slice_stats(self.current_stage)

    @property
    def filter_stats(self) -> pl.DataFrame:
        return (
            pl.DataFrame(stat.model_dump() for stat in self.filtering_stats).rename(
                {"n_fragments": "unique_fragments", "n_slices": "unique_slices"}
            )
            if self.filtering_stats
            else pl.DataFrame()
        )

    @property
    def read_stats(self) -> pl.DataFrame:
        filter_stats = self.filter_stats
        if filter_stats.is_empty():
            return pl.DataFrame()
        return (
            filter_stats.rename({"stage": "stat_type", "unique_fragments": "stat"})
            .select("stat_type", "stat")
            .with_columns(
                pl.lit("ccanalysis").alias("stage"),
                pl.lit(self.read_type).alias("read_type"),
                pl.lit(self.sample_name).alias("sample"),
                pl.lit(0).alias("read_number"),
            )
        )

    @property
    def captures(self) -> pl.DataFrame:
        return _captures(self.pipeline.slices)

    @property
    def slices_with_viewpoint(self) -> pl.DataFrame:
        return _slices_with_viewpoint(self.pipeline.slices)

    def filter_slices(
        self, output_slices: bool | str = False, output_location: Path | str = "."
    ) -> None:
        self.pipeline.run(output_slices=output_slices, output_location=output_location)
        self.filtering_stats = self.pipeline.filtering_stats
        self.current_stage = self.pipeline.current_stage
        self.current_filter = self.pipeline.current_filter

    def get_unfiltered_slices(self) -> None:
        self._apply(FilterStepName.GET_UNFILTERED_SLICES)

    def remove_unmapped_slices(self) -> None:
        self._apply(FilterStepName.REMOVE_UNMAPPED_SLICES)

    def remove_orphan_slices(self) -> None:
        self._apply(FilterStepName.REMOVE_ORPHAN_SLICES)

    def remove_duplicate_re_frags(self) -> None:
        self._apply(FilterStepName.REMOVE_DUPLICATE_RE_FRAGS)

    def remove_slices_without_re_frag_assigned(self) -> None:
        self._apply(FilterStepName.REMOVE_SLICES_WITHOUT_RE_FRAG_ASSIGNED)

    def remove_duplicate_slices(self) -> None:
        self._apply(FilterStepName.REMOVE_DUPLICATE_SLICES)

    def remove_duplicate_slices_pe(self) -> None:
        self._apply(FilterStepName.REMOVE_DUPLICATE_SLICES_PE)

    def remove_blacklisted_slices(self) -> None:
        self._apply(FilterStepName.REMOVE_BLACKLISTED_SLICES)

    def _apply(self, step_name: FilterStepName) -> None:
        self.pipeline.slices = FILTER_REGISTRY.get(step_name)(self.pipeline.slices)


class CCSliceFilter(SliceFilter):
    assay = Assay.CAPTURE

    @property
    def fragments(self) -> pl.DataFrame:
        return _capture_fragments(self.pipeline.slices)

    @property
    def reporters(self) -> pl.DataFrame:
        return self.pipeline.slices.filter(pl.col("capture_count") < 1)

    @property
    def capture_site_stats(self) -> pl.DataFrame:
        return self.captures["capture"].value_counts()

    @property
    def merged_captures_and_reporters(self) -> pl.DataFrame:
        captures = self.captures.rename(
            {
                column: f"capture_{column}"
                for column in self.captures.columns
                if column != "parent_read"
            }
        ).rename({"capture_capture": "capture"})
        reporters = self.reporters.rename(
            {
                column: f"reporter_{column}"
                for column in self.reporters.columns
                if column != "parent_read"
            }
        )
        return captures.join(reporters, on="parent_read", how="left")

    @property
    def cis_or_trans_stats(self) -> pl.DataFrame:
        return (
            self.merged_captures_and_reporters.with_columns(
                pl.when(pl.col("capture_chrom") == pl.col("reporter_chrom"))
                .then(pl.lit("cis"))
                .otherwise(pl.lit("trans"))
                .alias("cis/trans")
            )
            .group_by("capture", "cis/trans")
            .len(name="count")
            .rename({"capture": "viewpoint"})
            .with_columns(
                pl.lit(self.sample_name).alias("sample"),
                pl.lit(self.read_type).alias("read_type"),
            )
        )

    def remove_excluded_slices(self) -> None:
        self._apply(FilterStepName.REMOVE_EXCLUDED_SLICES)

    def remove_non_reporter_fragments(self) -> None:
        self._apply(FilterStepName.REMOVE_NON_REPORTER_FRAGMENTS)

    def remove_multi_capture_fragments(self) -> None:
        self._apply(FilterStepName.REMOVE_MULTI_CAPTURE_FRAGMENTS)

    def remove_viewpoint_adjacent_restriction_fragments(
        self, n_adjacent: int = 1
    ) -> None:
        self.pipeline.slices = remove_viewpoint_adjacent_restriction_fragments(
            self.pipeline.slices, n_adjacent=n_adjacent
        )


class TriCSliceFilter(CCSliceFilter):
    assay = Assay.TRI

    def remove_slices_with_one_reporter(self) -> None:
        self._apply(FilterStepName.REMOVE_SLICES_WITH_ONE_REPORTER)


class TiledCSliceFilter(SliceFilter):
    assay = Assay.TILED

    @property
    def captures(self) -> pl.DataFrame:
        return _captures(self.pipeline.slices)

    @property
    def fragments(self) -> pl.DataFrame:
        return _tiled_fragments(self.pipeline.slices)

    @property
    def cis_or_trans_stats(self) -> pl.DataFrame:
        rows = []
        slices = self.pipeline.slices
        capture_sites = (
            slices.filter(pl.col("capture_count") == 1)
            .select("capture")
            .drop_nulls()
            .unique(maintain_order=True)
            .get_column("capture")
            .to_list()
        )
        for capture_site in capture_sites:
            df_cap = slices.filter(
                (pl.col("capture_count") == 1) & (pl.col("capture") == capture_site)
            )
            if df_cap.is_empty():
                continue
            capture_chrom = df_cap.get_column("chrom").item(0)
            primary_slice_names = (
                df_cap.group_by("parent_read", maintain_order=True)
                .agg(pl.col("slice_name").first())
                .get_column("slice_name")
            )
            parent_reads = df_cap.select("parent_read").unique()
            df_not_primary_capture = df_cap.filter(
                ~pl.col("slice_name").is_in(primary_slice_names)
            )
            df_outside_capture = slices.filter(pl.col("capture_count") == 0).join(
                parent_reads, on="parent_read", how="semi"
            )
            df_pseudo_reporters = pl.concat(
                [df_not_primary_capture, df_outside_capture], how="vertical_relaxed"
            )
            n_cis_interactions = df_pseudo_reporters.filter(
                pl.col("chrom") == capture_chrom
            ).height
            n_trans_interactions = df_pseudo_reporters.height - n_cis_interactions
            rows.extend(
                [
                    {
                        "viewpoint": capture_site,
                        "cis/trans": "cis",
                        "count": n_cis_interactions,
                    },
                    {
                        "viewpoint": capture_site,
                        "cis/trans": "trans",
                        "count": n_trans_interactions,
                    },
                ]
            )

        if not rows:
            return pl.DataFrame(
                schema={
                    "viewpoint": pl.String,
                    "cis/trans": pl.String,
                    "count": pl.Int64,
                    "sample": pl.String,
                    "read_type": pl.String,
                }
            )

        return (
            pl.DataFrame(rows)
            .sort("viewpoint")
            .with_columns(
                pl.lit(self.sample_name).alias("sample"),
                pl.lit(self.read_type).alias("read_type"),
            )
        )

    def remove_slices_outside_capture(self) -> None:
        self._apply(FilterStepName.REMOVE_SLICES_OUTSIDE_CAPTURE)

    def remove_non_capture_fragments(self) -> None:
        self._apply(FilterStepName.REMOVE_NON_CAPTURE_FRAGMENTS)

    def remove_dual_capture_fragments(self) -> None:
        self._apply(FilterStepName.REMOVE_DUAL_CAPTURE_FRAGMENTS)

    def remove_religation(self) -> None:
        self._apply(FilterStepName.REMOVE_RELIGATION)
