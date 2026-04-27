from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import cast
import tomllib

from capcruncher.types import Assay, VALID_ASSAYS, validate_choice

from capcruncher.api.filtering.steps import (
    FilterStepInput,
    FilterStepName,
    _coerce_filter_step_name,
)

@dataclass(frozen=True)
class FilterStage:
    name: str
    steps: tuple[FilterStepName, ...]


@dataclass(frozen=True)
class FilterPlan:
    assay: Assay
    stages: tuple[FilterStage, ...]

    @classmethod
    def from_mapping(
        cls, data: Mapping[str, object], *, expected_assay: Assay | str | None = None
    ) -> "FilterPlan":
        assay_raw = data.get("assay")
        if not isinstance(assay_raw, str):
            raise ValueError("Filter profile must define string field 'assay'.")

        assay = validate_choice(assay_raw, VALID_ASSAYS, "assay")
        if expected_assay is not None:
            expected = validate_choice(expected_assay, VALID_ASSAYS, "expected_assay")
            if assay != expected:
                raise ValueError(
                    f"Filter profile assay {assay.value!r} does not match requested "
                    f"assay {expected.value!r}."
                )

        raw_stages = data.get("stages")
        if not isinstance(raw_stages, list) or not raw_stages:
            raise ValueError("Filter profile must define at least one [[stages]] entry.")

        stages: list[FilterStage] = []
        seen_stage_names: set[str] = set()
        for index, raw_stage in enumerate(raw_stages, start=1):
            if not isinstance(raw_stage, Mapping):
                raise ValueError(f"Filter stage {index} must be a table.")
            stage = cast(Mapping[str, object], raw_stage)

            name = stage.get("name")
            if not isinstance(name, str) or not name.strip():
                raise ValueError(f"Filter stage {index} must define a non-empty name.")
            if name in seen_stage_names:
                raise ValueError(f"Duplicate filter stage name: {name!r}.")
            seen_stage_names.add(name)

            raw_steps = stage.get("steps")
            if not isinstance(raw_steps, list) or not raw_steps:
                raise ValueError(f"Filter stage {name!r} must define non-empty steps.")
            if not all(isinstance(step, str) and step.strip() for step in raw_steps):
                raise ValueError(f"Filter stage {name!r} contains an invalid step name.")
            step_names = cast(list[str], raw_steps)

            stages.append(
                FilterStage(
                    name=name,
                    steps=tuple(_coerce_filter_step_name(step) for step in step_names),
                )
            )

        return cls(assay=assay, stages=tuple(stages))

    @classmethod
    def from_toml(
        cls, path: Path | str, *, expected_assay: Assay | str | None = None
    ) -> "FilterPlan":
        profile_path = Path(path)
        if profile_path.suffix.lower() in {".yaml", ".yml"}:
            raise ValueError(
                "YAML filter profiles are no longer supported. Use a TOML filter "
                "profile with [[stages]] entries."
            )
        with profile_path.open("rb") as handle:
            data = tomllib.load(handle)
        return cls.from_mapping(data, expected_assay=expected_assay)



DEFAULT_FILTER_PLANS: dict[Assay, FilterPlan] = {
    Assay.CAPTURE: FilterPlan(
        assay=Assay.CAPTURE,
        stages=(
            FilterStage("pre-filtering", (FilterStepName.GET_UNFILTERED_SLICES,)),
            FilterStage("mapped", (FilterStepName.REMOVE_UNMAPPED_SLICES,)),
            FilterStage(
                "contains_single_capture",
                (
                    FilterStepName.REMOVE_ORPHAN_SLICES,
                    FilterStepName.REMOVE_MULTI_CAPTURE_FRAGMENTS,
                ),
            ),
            FilterStage(
                "contains_capture_and_reporter",
                (
                    FilterStepName.REMOVE_EXCLUDED_SLICES,
                    FilterStepName.REMOVE_BLACKLISTED_SLICES,
                    FilterStepName.REMOVE_NON_REPORTER_FRAGMENTS,
                    FilterStepName.REMOVE_VIEWPOINT_ADJACENT_RESTRICTION_FRAGMENTS,
                ),
            ),
            FilterStage(
                "duplicate_filtered",
                (
                    FilterStepName.REMOVE_SLICES_WITHOUT_RE_FRAG_ASSIGNED,
                    FilterStepName.REMOVE_DUPLICATE_RE_FRAGS,
                    FilterStepName.REMOVE_DUPLICATE_SLICES,
                    FilterStepName.REMOVE_DUPLICATE_SLICES_PE,
                    FilterStepName.REMOVE_NON_REPORTER_FRAGMENTS,
                ),
            ),
        ),
    ),
    Assay.TRI: FilterPlan(
        assay=Assay.TRI,
        stages=(
            FilterStage("pre-filtering", (FilterStepName.GET_UNFILTERED_SLICES,)),
            FilterStage(
                "mapped",
                (
                    FilterStepName.REMOVE_UNMAPPED_SLICES,
                    FilterStepName.REMOVE_SLICES_WITHOUT_RE_FRAG_ASSIGNED,
                ),
            ),
            FilterStage(
                "contains_single_capture",
                (
                    FilterStepName.REMOVE_ORPHAN_SLICES,
                    FilterStepName.REMOVE_MULTI_CAPTURE_FRAGMENTS,
                ),
            ),
            FilterStage(
                "contains_capture_and_reporter",
                (
                    FilterStepName.REMOVE_BLACKLISTED_SLICES,
                    FilterStepName.REMOVE_NON_REPORTER_FRAGMENTS,
                ),
            ),
            FilterStage(
                "duplicate_filtered",
                (
                    FilterStepName.REMOVE_DUPLICATE_RE_FRAGS,
                    FilterStepName.REMOVE_DUPLICATE_SLICES,
                    FilterStepName.REMOVE_DUPLICATE_SLICES_PE,
                    FilterStepName.REMOVE_NON_REPORTER_FRAGMENTS,
                ),
            ),
            FilterStage(
                "tric_reporter", (FilterStepName.REMOVE_SLICES_WITH_ONE_REPORTER,)
            ),
        ),
    ),
    Assay.TILED: FilterPlan(
        assay=Assay.TILED,
        stages=(
            FilterStage("pre-filtering", (FilterStepName.GET_UNFILTERED_SLICES,)),
            FilterStage(
                "mapped",
                (
                    FilterStepName.REMOVE_UNMAPPED_SLICES,
                    FilterStepName.REMOVE_ORPHAN_SLICES,
                ),
            ),
            FilterStage("not_blacklisted", (FilterStepName.REMOVE_BLACKLISTED_SLICES,)),
            FilterStage(
                "contains_capture",
                (
                    FilterStepName.REMOVE_NON_CAPTURE_FRAGMENTS,
                    FilterStepName.REMOVE_DUAL_CAPTURE_FRAGMENTS,
                ),
            ),
            FilterStage(
                "duplicate_filtered",
                (
                    FilterStepName.REMOVE_SLICES_WITHOUT_RE_FRAG_ASSIGNED,
                    FilterStepName.REMOVE_DUPLICATE_RE_FRAGS,
                    FilterStepName.REMOVE_DUPLICATE_SLICES,
                    FilterStepName.REMOVE_DUPLICATE_SLICES_PE,
                ),
            ),
            FilterStage(
                "has_reporter",
                (FilterStepName.REMOVE_ORPHAN_SLICES, FilterStepName.REMOVE_RELIGATION),
            ),
        ),
    ),
}


