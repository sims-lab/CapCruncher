import sys
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pyranges1 as pr
from loguru import logger
from pydantic import BaseModel, PositiveFloat, PositiveInt, field_validator, model_validator

from capcruncher.api.intervals.annotate import annotate_intervals, remove_duplicates_from_bed
from capcruncher.types import (
    AnnotationAction,
    DuplicateAction,
    InvalidBedAction,
    VALID_ANNOTATION_ACTIONS,
    VALID_DUPLICATE_ACTIONS,
    VALID_INVALID_BED_ACTIONS,
    validate_choice,
    validate_choices,
)
from capcruncher.utils import (
    convert_bed_to_pr,
    cycle_argument,
    hash_column,
)

warnings.simplefilter("ignore")


def _bam_to_bed_dataframe(bam_path: Path | str) -> pd.DataFrame:
    import pysam

    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_name is None:
                continue

            rows.append(
                {
                    "chrom": read.reference_name,
                    "start": read.reference_start,
                    "end": read.reference_end,
                    "name": read.query_name,
                    "score": read.mapping_quality,
                    "strand": "-" if read.is_reverse else "+",
                }
            )

    return pd.DataFrame.from_records(
        rows, columns=["chrom", "start", "end", "name", "score", "strand"]
    )


class AlignmentAnnotateOptions(BaseModel):
    """Validated options for alignment annotation."""

    slices: Path | str
    actions: Sequence[AnnotationAction | str] | None = ()
    bed_files: Sequence[Path | str] | None = ()
    names: Sequence[str] | None = ()
    overlap_fractions: Sequence[PositiveFloat | float] | None = (1e-9,)
    output: Path = Path("annotated.slices.parquet")
    duplicates: DuplicateAction | str = DuplicateAction.REMOVE
    invalid_bed_action: InvalidBedAction | str = InvalidBedAction.ERROR
    n_cores: PositiveInt = 1
    blacklist: Path | None = None
    prioritize_cis_slices: bool = False
    priority_chroms: Sequence[str] | str | None = ()

    @field_validator("actions", mode="before")
    @classmethod
    def validate_actions(
        cls, value: Sequence[str | AnnotationAction] | None
    ) -> tuple[AnnotationAction, ...]:
        return validate_choices(tuple(value or ()), VALID_ANNOTATION_ACTIONS, "actions")

    @field_validator("bed_files", mode="before")
    @classmethod
    def validate_bed_files(cls, value: Sequence[Path | str] | None) -> tuple[Path, ...]:
        return tuple(Path(path) for path in (value or ()))

    @field_validator("names", mode="before")
    @classmethod
    def validate_names(cls, value: Sequence[str] | None) -> tuple[str, ...]:
        return tuple(value or ())

    @field_validator("overlap_fractions", mode="before")
    @classmethod
    def validate_overlap_fractions(
        cls, value: Sequence[float] | None
    ) -> tuple[float, ...]:
        return tuple(value or (1e-9,))

    @field_validator("slices", mode="before")
    @classmethod
    def validate_slices(cls, value: Path | str) -> Path | str:
        if value == "-":
            return value
        path = Path(value)
        if not path.exists():
            raise ValueError(f"slices does not exist: {path}")
        return path

    @field_validator("output")
    @classmethod
    def validate_output_parent(cls, value: Path) -> Path:
        if value.parent and not value.parent.exists():
            raise ValueError(f"output parent directory does not exist: {value.parent}")
        return value

    @field_validator("duplicates", mode="before")
    @classmethod
    def validate_duplicates(cls, value: str | DuplicateAction) -> DuplicateAction:
        return validate_choice(value, VALID_DUPLICATE_ACTIONS, "duplicates")

    @field_validator("invalid_bed_action", mode="before")
    @classmethod
    def validate_invalid_bed_action(
        cls, value: str | InvalidBedAction
    ) -> InvalidBedAction:
        return validate_choice(value, VALID_INVALID_BED_ACTIONS, "invalid_bed_action")

    @field_validator("blacklist", mode="before")
    @classmethod
    def validate_blacklist(cls, value: Path | str | None) -> Path | None:
        if value in (None, ""):
            return None
        return Path(value)

    @field_validator("priority_chroms", mode="before")
    @classmethod
    def validate_priority_chroms(
        cls, value: Sequence[str] | str | None
    ) -> tuple[str, ...]:
        if value in (None, ""):
            return ()
        if isinstance(value, str):
            return tuple(chrom for chrom in value.split(",") if chrom)
        return tuple(value)

    @model_validator(mode="after")
    def validate_annotation_lengths(self) -> "AlignmentAnnotateOptions":
        actions = tuple(self.actions or ())
        bed_files = tuple(self.bed_files or ())
        names = tuple(self.names or ())
        if len(actions) != len(bed_files) or len(names) != len(bed_files):
            raise ValueError(
                "The lengths of the supplied bed files, actions, and names do not match."
            )
        self.actions = actions
        self.bed_files = bed_files
        self.names = names
        self.overlap_fractions = tuple(self.overlap_fractions or (1e-9,))
        self.priority_chroms = tuple(self.priority_chroms or ())
        return self


def annotate(
    slices: Path | str,
    actions: Sequence[AnnotationAction | str] | None = None,
    bed_files: Sequence[Path | str] | None = None,
    names: Sequence[str] | None = None,
    overlap_fractions: Sequence[float] | None = None,
    output: Path | str | None = None,
    duplicates: DuplicateAction | str = DuplicateAction.REMOVE,
    n_cores: int = 1,
    blacklist: Path | str | None = None,
    prioritize_cis_slices: bool = False,
    priority_chroms: Sequence[str] | str | None = None,
    invalid_bed_action: InvalidBedAction | str = InvalidBedAction.ERROR,
    **kwargs: Any,
) -> None:
    """Annotate a BED-like input with one or more BED files.

    Args:
        slices: Input BED path, BAM path, or ``"-"`` to read BED rows from stdin.
        actions: Annotation methods. Valid values are ``"get"`` and ``"count"``.
        bed_files: BED files to intersect with ``slices``.
        names: Output annotation column names. Must match ``bed_files`` length.
        overlap_fractions: Minimum overlap fractions used for intersections.
        output: Output parquet path.
        duplicates: Duplicate reconciliation method. Only ``"remove"`` is supported.
        n_cores: Number of cores requested for annotation.
        blacklist: Optional BED file of regions to subtract before annotation.
        prioritize_cis_slices: Prefer cis slices while removing duplicates.
        priority_chroms: Chromosomes to prioritize during duplicate removal.
        invalid_bed_action: How invalid annotation BED files are handled.

    Raises:
        ValueError: If constrained string options or paired option lengths are invalid.
    """

    options = AlignmentAnnotateOptions(
        slices=slices,
        actions=actions,
        bed_files=bed_files,
        names=names,
        overlap_fractions=overlap_fractions,
        output=Path(output) if output is not None else Path("annotated.slices.parquet"),
        duplicates=duplicates,
        n_cores=n_cores,
        blacklist=Path(blacklist) if blacklist not in (None, "") else None,
        prioritize_cis_slices=prioritize_cis_slices,
        priority_chroms=priority_chroms,
        invalid_bed_action=invalid_bed_action,
    )

    with logger.catch():
        logger.info("Validating commandline arguments")
        slices_input = options.slices

        if slices_input == "-":
            logger.info("Reading slices from stdin")
            slice_intervals = convert_bed_to_pr(
                pd.read_csv(sys.stdin, sep="\t", header=None)
            )

        elif str(slices_input).endswith(".bam"):
            logger.info("Converting bam to bed")
            slice_intervals = _bam_to_bed_dataframe(slices_input).pipe(
                convert_bed_to_pr
            )

        else:
            slice_intervals = convert_bed_to_pr(slices_input)

        logger.info("Validating input bed file before annotation")

        if options.blacklist:
            try:
                logger.info("Removing blacklisted regions from the bed file")
                gr_blacklist = pr.read_bed(options.blacklist)
                slice_intervals = slice_intervals.subtract_overlaps(
                    gr_blacklist, strand_behavior="ignore"
                )
            except Exception as e:
                logger.warning(
                    f"Failed to remove blacklisted regions from the bed file. {e}"
                )

        logger.info("Dealing with duplicates in the bed file")

        slice_intervals = remove_duplicates_from_bed(
            slice_intervals,
            prioritize_cis_slices=options.prioritize_cis_slices,
            chroms_to_prioritize=list(options.priority_chroms or ()) or None,
        )

        actions = cast(tuple[AnnotationAction, ...], options.actions or ())
        bed_files = tuple(options.bed_files or ())
        names = tuple(options.names or ())
        overlap_fractions = tuple(options.overlap_fractions or (1e-9,))
        for action, bed_file, name, fraction in zip(
            actions,
            bed_files,
            names,
            cycle_argument(overlap_fractions),
        ):
            logger.info(
                f"Performing {name} intersection with {bed_file} using {action} method with {fraction} overlap fraction. {len(slice_intervals)} slices to intersect."
            )

            slice_intervals = annotate_intervals(
                query=slice_intervals,
                annotations=bed_file,
                name=name,
                method=action,
                fraction=fraction,
            )

        logger.info("Writing annotations to file.")
        df_annotation = slice_intervals.rename(columns={"Name": "slice_name"}).assign(
            slice_id=lambda df: hash_column(df.slice_name)
        )
        df_annotation.to_parquet(options.output, compression="snappy")
