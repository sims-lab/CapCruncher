from __future__ import annotations

import re
from typing import Annotated, Any

from pydantic import BaseModel, BeforeValidator, ConfigDict, Field, field_validator, model_validator

from capcruncher.types import (
    Assay,
    FastqSplitMethod,
    FLAG_OFF_VALUES,
    FLAG_ON_VALUES,
    FLAG_NONE_VALUES,
    SummaryMethod,
    VALID_ASSAYS,
    VALID_FASTQ_SPLIT_METHODS,
    VALID_SUMMARY_METHODS,
    validate_choice,
)
from capcruncher.utils import get_restriction_site


HUB_COLOR_BY_ALIASES = {
    "samplename": "sample",
    "sample_name": "sample",
    "method": "aggregation",
}
HUB_COLOR_BY_COLUMNS = {
    "aggregation",
    "category",
    "normalisation",
    "overlay",
    "sample",
    "viewpoint",
}

PLOT_NORMALISATIONS = {
    "raw",
    "n_interactions",
    "n_rf_n_interactions",
    "ice",
    "icen_cis",
    "icen_scale",
}
ALIGNERS = {"bowtie", "bowtie2"}
PRIORITY_CHROMOSOME_MODES = {"viewpoints"}


def _coerce_flag_value(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    s = str(v).strip().lower()
    if s in FLAG_ON_VALUES:
        return True
    if s in FLAG_OFF_VALUES | FLAG_NONE_VALUES:
        return False
    return v


FlagValue = Annotated[bool, BeforeValidator(_coerce_flag_value)]


def is_on(param: Any) -> bool:
    return str(param).strip().lower() in FLAG_ON_VALUES


def is_off(param: Any) -> bool:
    return str(param).strip().lower() in FLAG_OFF_VALUES


def is_none(param: Any) -> bool:
    return str(param).strip().lower() in FLAG_NONE_VALUES


def normalise_scalar_config_value(value: Any) -> Any:
    if isinstance(value, bool):
        return value

    if value is None:
        return False

    if not isinstance(value, str):
        return value

    value = value.strip()
    if is_on(value):
        return True
    if is_off(value) or is_none(value):
        return False
    return value


def normalise_config_values(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: normalise_config_values(entry) for key, entry in value.items()}

    if isinstance(value, list):
        return [normalise_config_values(entry) for entry in value]

    return normalise_scalar_config_value(value)


def normalise_hub_color_by(value: str | bool | None) -> str | None:
    """Map hub color fields onto known TrackNado metadata columns."""
    if value is None or value is False:
        return None

    value = str(value).strip()
    if value.lower() in {"", "none", "false", "no", "null"}:
        return None

    return HUB_COLOR_BY_ALIASES.get(value.lower(), value)


class PipelineBaseModel(BaseModel):
    model_config = ConfigDict(extra="allow", use_enum_values=True)


class AnalysisConfig(PipelineBaseModel):
    method: Assay | str | None = None
    viewpoints: str | bool | None = None
    restriction_enzyme: str | None = None
    bin_sizes: list[int] = Field(default_factory=list)
    reporter_exclusion_zone: int | None = None
    design: str | bool | None = None
    regenerate_fastq: bool | None = None
    blacklist: str | bool | None = None
    filter_profile: str | bool | None = None

    @field_validator("method", mode="before")
    @classmethod
    def validate_method(cls, value: str | bool | None) -> str | None:
        if value is None or value is False:
            return None
        return validate_choice(str(value).lower(), VALID_ASSAYS, "analysis.method").value

    @field_validator("bin_sizes", mode="before")
    @classmethod
    def validate_bin_sizes(cls, value: Any) -> list[int]:
        if value is None or value is False:
            return []

        if isinstance(value, int):
            values = [value]
        elif isinstance(value, str):
            values = [entry for entry in re.split(r"[,;]\s*|\s+", value) if entry]
        elif isinstance(value, list):
            values = value
        else:
            raise ValueError(
                "analysis.bin_sizes must be an int, list of ints, or separated string."
            )

        bins = [int(entry) for entry in values]
        if any(bin_size < 0 for bin_size in bins):
            raise ValueError("analysis.bin_sizes cannot contain negative values.")
        return bins

    @field_validator("restriction_enzyme")
    @classmethod
    def validate_restriction_enzyme(cls, value: str | None) -> str | None:
        if value is None:
            return value

        try:
            get_restriction_site(value)
        except ValueError as exc:
            raise ValueError(
                "analysis.restriction_enzyme must be a known enzyme name or an "
                f"explicit DNA recognition sequence. Got: {value!r}"
            ) from exc
        return value

    @field_validator("reporter_exclusion_zone")
    @classmethod
    def validate_reporter_exclusion_zone(cls, value: int | None) -> int | None:
        if value is not None and value < 0:
            raise ValueError("analysis.reporter_exclusion_zone cannot be negative.")
        return value


class GenomeConfig(PipelineBaseModel):
    name: str | None = None
    fasta: str | bool | None = None
    aligner_index: str | bool | None = None
    chrom_sizes: str | bool | None = None
    custom: bool | None = None
    organism: str | bool | None = None
    twobit: str | bool | None = None
    genome_default_position: str | bool | None = None


class HubConfig(PipelineBaseModel):
    create: bool | None = None
    dir: str | bool | None = None
    name: str | bool | None = None
    email: str | bool | None = None
    color_by: str | bool | None = "sample"
    custom_genome: bool | None = None

    @field_validator("color_by", mode="before")
    @classmethod
    def validate_color_by(cls, value: str | bool | None) -> str | bool:
        color_by = normalise_hub_color_by(value)
        if color_by is None:
            return False

        if color_by in HUB_COLOR_BY_COLUMNS:
            return color_by

        valid_values = sorted(
            HUB_COLOR_BY_COLUMNS | set(HUB_COLOR_BY_ALIASES) | {"none"}
        )
        raise ValueError(
            f"Invalid hub.color_by value {value!r}. "
            f"Choose one of: {', '.join(valid_values)}."
        )


class PlotConfig(PipelineBaseModel):
    create: bool | None = None
    coordinates: str | bool | None = None
    normalisation: str | bool | None = None
    genes: str | bool | None = None

    @field_validator("normalisation", mode="before")
    @classmethod
    def validate_normalisation(cls, value: str | bool | None) -> str | bool | None:
        if value is None or value is False:
            return value

        value = str(value)
        if value not in PLOT_NORMALISATIONS:
            raise ValueError(
                "plot.normalisation must be one of: "
                f"{', '.join(sorted(PLOT_NORMALISATIONS))}. Got: {value!r}"
            )
        return value


class AlignConfig(PipelineBaseModel):
    aligner: str | None = None
    index_flag: str | bool | None = None
    options: str | bool | None = None

    @field_validator("aligner")
    @classmethod
    def validate_aligner(cls, value: str | None) -> str | None:
        if value is None:
            return value

        value = value.lower()
        if value not in ALIGNERS:
            raise ValueError(
                f"align.aligner must be one of: {', '.join(sorted(ALIGNERS))}. "
                f"Got: {value!r}"
            )
        return value


class AnalysisOptionalConfig(PipelineBaseModel):
    filter_profile: str | bool | None = None
    blacklist: str | bool | None = None
    prioritize_cis_slices: bool | None = None
    priority_chromosomes: str | bool | None = None
    minimum_viewpoint_overlap: float | None = None
    force_bigwig_generation: bool | None = None

    @field_validator("priority_chromosomes", mode="before")
    @classmethod
    def validate_priority_chromosomes(cls, value: str | bool | None) -> str | bool:
        if value is None or value is False:
            return False

        value = str(value).strip()
        chromosomes = [chromosome.strip() for chromosome in value.split(",")]
        if value in PRIORITY_CHROMOSOME_MODES or all(chromosomes):
            return value

        raise ValueError(
            "analysis_optional.priority_chromosomes must be 'viewpoints', "
            "a comma-separated chromosome list, or none."
        )

    @field_validator("minimum_viewpoint_overlap")
    @classmethod
    def validate_minimum_viewpoint_overlap(
        cls, value: float | None
    ) -> float | None:
        if value is not None and not 0 <= value <= 1:
            raise ValueError(
                "analysis_optional.minimum_viewpoint_overlap must be between 0 and 1."
            )
        return value


class NormalisationConfig(PipelineBaseModel):
    scale_factor: int | None = None
    regions: str | bool | None = None

    @field_validator("scale_factor")
    @classmethod
    def validate_scale_factor(cls, value: int | None) -> int | None:
        if value is not None and value <= 0:
            raise ValueError("normalisation.scale_factor must be greater than 0.")
        return value


class DifferentialConfig(PipelineBaseModel):
    contrast: str | bool | None = None
    distance: int | None = None

    @field_validator("distance")
    @classmethod
    def validate_distance(cls, value: int | None) -> int | None:
        if value is not None and value < 0:
            raise ValueError("differential.distance cannot be negative.")
        return value


class SplitConfig(PipelineBaseModel):
    n_reads: int | None = None
    method: FastqSplitMethod | str | None = None

    @field_validator("n_reads")
    @classmethod
    def validate_n_reads(cls, value: int | None) -> int | None:
        if value is not None and value <= 0:
            raise ValueError("split.n_reads must be greater than 0.")
        return value

    @field_validator("method", mode="before")
    @classmethod
    def validate_method(cls, value: str | bool | None) -> str | None:
        if value is None or value is False:
            return None
        return validate_choice(
            str(value).lower(), VALID_FASTQ_SPLIT_METHODS, "split.method"
        ).value


class CompareConfig(PipelineBaseModel):
    summary_methods: str | None = None

    @field_validator("summary_methods", mode="before")
    @classmethod
    def validate_summary_methods(cls, value: str | bool | None) -> str | None:
        if value is None or value is False:
            return None

        values = [entry for entry in re.split(r"[,;]\s*|\s+", str(value)) if entry]
        methods = [
            validate_choice(method, VALID_SUMMARY_METHODS, "compare.summary_methods")
            for method in values
        ]
        return ",".join(method.value for method in methods)


class ExecutionConfig(PipelineBaseModel):
    container_image: str | bool | None = None


class TrimConfig(PipelineBaseModel):
    options: str | bool | None = None


class PipelineConfig(PipelineBaseModel):
    """Permissive model for normalising workflow config before Snakemake use."""

    version: str | None = None
    analysis: AnalysisConfig | None = None
    analysis_optional: AnalysisOptionalConfig | None = None
    align: AlignConfig | None = None
    compare: CompareConfig | None = None
    execution: ExecutionConfig | None = None
    genome: GenomeConfig | None = None
    hub: HubConfig | None = None
    plot: PlotConfig | None = None
    normalisation: NormalisationConfig | None = None
    differential: DifferentialConfig | None = None
    split: SplitConfig | None = None
    trim: TrimConfig | None = None

    @model_validator(mode="after")
    def validate_custom_hub_genome(self) -> "PipelineConfig":
        if not self.hub:
            return self

        custom_genome = bool(self.hub.custom_genome)
        if self.genome:
            custom_genome = custom_genome or bool(self.genome.custom)

        if not (self.hub.create and custom_genome):
            return self

        if not (self.genome and self.genome.twobit):
            raise ValueError(
                "Custom UCSC hub genomes require genome.twobit when hub.create is true."
            )

        self.hub.custom_genome = True
        return self


def validate_pipeline_config(config: dict[str, Any]) -> dict[str, Any]:
    formatted = normalise_config_values(config)
    validated = PipelineConfig.model_validate(formatted)
    return validated.model_dump(mode="python", exclude_none=True)


def format_pipeline_config(config: dict[str, Any]) -> dict[str, Any]:
    """Normalise and validate the pipeline config in place."""
    validated = validate_pipeline_config(config)
    config.clear()
    config.update(validated)
    return config
