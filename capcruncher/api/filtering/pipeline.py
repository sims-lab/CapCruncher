from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import ClassVar

import numpy as np
import pandas as pd
import polars as pl
from loguru import logger

from capcruncher.api.statistics import SliceFilterStats
from capcruncher.types import Assay, ReadType

from capcruncher.api.filtering.plan import DEFAULT_FILTER_PLANS, FilterPlan, FilterStage
from capcruncher.api.filtering.steps import (
    FILTER_REGISTRY,
    FilterStepInput,
    FilterStepName,
    FilterStepRegistry,
    _captures,
    _normalise_slices,
    _slices_with_viewpoint,
    _capture_fragments,
    _tiled_fragments,
    remove_viewpoint_adjacent_restriction_fragments,
)

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
        filter_stages: FilterPlan | Mapping[str, Sequence[FilterStepInput]] | None = None,
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
        self._filter_stats = pd.DataFrame()

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
    def slices(self) -> pd.DataFrame:
        return self.pipeline.slices.to_pandas()

    @slices.setter
    def slices(self, slices: pd.DataFrame | pl.DataFrame) -> None:
        self.pipeline.slices = _normalise_slices(slices)

    @property
    def filters(self) -> list[str]:
        return FILTER_REGISTRY.names

    @property
    def filter_stages(self) -> dict[str, list[str]]:
        return {stage.name: [step.value for step in stage.steps] for stage in self.plan.stages}

    @property
    def slice_stats(self) -> SliceFilterStats:
        return self.pipeline.slice_stats(self.current_stage)

    @property
    def filter_stats(self) -> pd.DataFrame:
        return (
            pd.DataFrame(stat.model_dump() for stat in self.filtering_stats)
            .rename(columns={"n_fragments": "unique_fragments", "n_slices": "unique_slices"})
            if self.filtering_stats
            else pd.DataFrame()
        )

    @property
    def read_stats(self) -> pd.DataFrame:
        if self.filter_stats.empty:
            return pd.DataFrame()
        return self.filter_stats.rename(
            columns={"stage": "stat_type", "unique_fragments": "stat"}
        )[["stat_type", "stat"]].assign(
            stage="ccanalysis",
            read_type=self.read_type,
            sample=self.sample_name,
            read_number=0,
        )

    @property
    def captures(self) -> pd.DataFrame:
        return _captures(self.pipeline.slices).to_pandas()

    @property
    def slices_with_viewpoint(self) -> pd.DataFrame:
        return _slices_with_viewpoint(self.pipeline.slices).to_pandas()

    def filter_slices(self, output_slices: bool | str = False, output_location: Path | str = ".") -> None:
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
    def fragments(self) -> pd.DataFrame:
        return _capture_fragments(self.pipeline.slices).to_pandas()

    @property
    def reporters(self) -> pd.DataFrame:
        return self.pipeline.slices.filter(pl.col("capture_count") < 1).to_pandas()

    @property
    def capture_site_stats(self) -> pd.Series:
        return self.captures["capture"].value_counts()

    @property
    def merged_captures_and_reporters(self) -> pd.DataFrame:
        captures = (
            self.captures.set_index("parent_read")
            .add_prefix("capture_")
            .rename(columns={"capture_capture": "capture"})
        )
        reporters = self.reporters.set_index("parent_read").add_prefix("reporter_")
        return captures.join(reporters).reset_index()

    @property
    def cis_or_trans_stats(self) -> pd.DataFrame:
        cap_and_rep = self.merged_captures_and_reporters.copy()
        cap_and_rep["cis/trans"] = np.where(
            cap_and_rep["capture_chrom"] == cap_and_rep["reporter_chrom"],
            "cis",
            "trans",
        )
        return (
            cap_and_rep.groupby(["capture", "cis/trans"])
            .size()
            .reset_index()
            .rename(columns={"capture": "viewpoint", 0: "count"})
            .assign(sample=self.sample_name, read_type=self.read_type)
        )

    def remove_excluded_slices(self) -> None:
        self._apply(FilterStepName.REMOVE_EXCLUDED_SLICES)

    def remove_non_reporter_fragments(self) -> None:
        self._apply(FilterStepName.REMOVE_NON_REPORTER_FRAGMENTS)

    def remove_multi_capture_fragments(self) -> None:
        self._apply(FilterStepName.REMOVE_MULTI_CAPTURE_FRAGMENTS)

    def remove_viewpoint_adjacent_restriction_fragments(self, n_adjacent: int = 1) -> None:
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
    def captures(self) -> pd.DataFrame:
        return _captures(self.pipeline.slices).to_pandas()

    @property
    def fragments(self) -> pd.DataFrame:
        return _tiled_fragments(self.pipeline.slices).to_pandas()

    @property
    def cis_or_trans_stats(self) -> pd.DataFrame:
        interactions_by_capture = {}
        for capture_site, df_cap in self.slices.query("capture_count == 1").groupby(
            "capture", observed=False
        ):
            capture_chrom = df_cap.iloc[0]["chrom"]
            df_primary_capture = df_cap.groupby("parent_read").first()
            df_not_primary_capture = df_cap.loc[
                ~(df_cap["slice_name"].isin(df_primary_capture["slice_name"]))
            ]
            df_outside_capture = self.slices.query("capture_count == 0").loc[
                lambda df_rep: df_rep["parent_read"].isin(df_cap["parent_read"])
            ]
            df_pseudo_reporters = pd.concat(
                [df_not_primary_capture, df_outside_capture]
            )
            n_cis_interactions = df_pseudo_reporters.query(
                f'chrom == "{capture_chrom}"'
            ).shape[0]
            n_trans_interactions = df_pseudo_reporters.shape[0] - n_cis_interactions
            interactions_by_capture[capture_site] = {
                "cis": n_cis_interactions,
                "trans": n_trans_interactions,
            }

        return (
            pd.DataFrame(interactions_by_capture)
            .transpose()
            .reset_index()
            .rename(columns={"index": "capture"})
            .melt(id_vars="capture", var_name="cis/trans", value_name="count")
            .sort_values("capture")
            .assign(sample=self.sample_name, read_type=self.read_type)
            .rename(columns={"capture": "viewpoint"})
        )

    def remove_slices_outside_capture(self) -> None:
        self._apply(FilterStepName.REMOVE_SLICES_OUTSIDE_CAPTURE)

    def remove_non_capture_fragments(self) -> None:
        self._apply(FilterStepName.REMOVE_NON_CAPTURE_FRAGMENTS)

    def remove_dual_capture_fragments(self) -> None:
        self._apply(FilterStepName.REMOVE_DUAL_CAPTURE_FRAGMENTS)

    def remove_religation(self) -> None:
        self._apply(FilterStepName.REMOVE_RELIGATION)
