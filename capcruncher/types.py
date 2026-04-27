from collections.abc import Iterable
from enum import StrEnum
from pathlib import Path
from typing import TypeVar


class Assay(StrEnum):
    CAPTURE = "capture"
    TRI = "tri"
    TILED = "tiled"


class ReadType(StrEnum):
    FLASHED = "flashed"
    PE = "pe"


class AnnotationAction(StrEnum):
    GET = "get"
    COUNT = "count"


class DuplicateAction(StrEnum):
    REMOVE = "remove"


class InvalidBedAction(StrEnum):
    ERROR = "error"
    IGNORE = "ignore"


class Normalisation(StrEnum):
    RAW = "raw"
    N_CIS = "n_cis"
    REGION = "region"


class BedgraphFormat(StrEnum):
    BEDGRAPH = "bedgraph"
    BIGWIG = "bigwig"


class CompareFormat(StrEnum):
    AUTO = "auto"
    COOLER = "cooler"
    BEDGRAPH = "bedgraph"


class OutputFormat(StrEnum):
    BEDGRAPH = "bedgraph"
    TSV = "tsv"


class SummaryMethod(StrEnum):
    MEAN = "mean"


class DictFormat(StrEnum):
    JSON = "json"
    PICKLE = "pickle"


class DictDType(StrEnum):
    INT = "int"
    STR = "str"


class Executor(StrEnum):
    LOCAL = "local"
    PROCESS = "process"
    RAY = "ray"


class FastqSplitMethod(StrEnum):
    PYTHON = "python"
    UNIX = "unix"
    SEQKIT = "seqkit"


class FastqSplitType(StrEnum):
    N_READS = "n-reads"
    N_PARTS = "n-parts"


class CisOrTrans(StrEnum):
    CIS = "cis"
    TRANS = "trans"


class BinningMethod(StrEnum):
    OVERLAP = "overlap"
    MIDPOINT = "midpoint"


VALID_ASSAYS: tuple[Assay, ...] = tuple(Assay)
VALID_READ_TYPES: tuple[ReadType, ...] = tuple(ReadType)
VALID_ANNOTATION_ACTIONS: tuple[AnnotationAction, ...] = tuple(AnnotationAction)
VALID_DUPLICATE_ACTIONS: tuple[DuplicateAction, ...] = tuple(DuplicateAction)
VALID_INVALID_BED_ACTIONS: tuple[InvalidBedAction, ...] = tuple(InvalidBedAction)
VALID_NORMALISATIONS: tuple[Normalisation, ...] = tuple(Normalisation)
VALID_BEDGRAPH_FORMATS: tuple[BedgraphFormat, ...] = tuple(BedgraphFormat)
VALID_COMPARE_FORMATS: tuple[CompareFormat, ...] = tuple(CompareFormat)
VALID_OUTPUT_FORMATS: tuple[OutputFormat, ...] = tuple(OutputFormat)
VALID_SUMMARY_METHODS: tuple[SummaryMethod, ...] = tuple(SummaryMethod)
VALID_DICT_FORMATS: tuple[DictFormat, ...] = tuple(DictFormat)
VALID_DICT_DTYPES: tuple[DictDType, ...] = tuple(DictDType)
VALID_EXECUTORS: tuple[Executor, ...] = tuple(Executor)
VALID_FASTQ_SPLIT_METHODS: tuple[FastqSplitMethod, ...] = tuple(FastqSplitMethod)
VALID_FASTQ_SPLIT_TYPES: tuple[FastqSplitType, ...] = tuple(FastqSplitType)
VALID_CIS_OR_TRANS: tuple[CisOrTrans, ...] = tuple(CisOrTrans)
VALID_BINNING_METHODS: tuple[BinningMethod, ...] = tuple(BinningMethod)

Choice = TypeVar("Choice", bound=StrEnum)


def _choice_value(choice: str | StrEnum) -> str:
    return choice.value if isinstance(choice, StrEnum) else choice


def _format_choices(valid: Iterable[str | StrEnum]) -> str:
    return ", ".join(_choice_value(choice) for choice in valid)


def validate_choice(
    value: str | Choice, valid: tuple[Choice, ...], option_name: str
) -> Choice:
    """Return a string enum option value or raise a clear validation error."""
    if isinstance(value, StrEnum) and value in valid:
        return value

    value_str = str(value)
    try:
        enum_type = type(valid[0])
        return enum_type(value_str)
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"{option_name} must be one of: {_format_choices(valid)}. Got: {value_str!r}"
        ) from exc


def validate_choices(
    values: Iterable[str | Choice], valid: tuple[Choice, ...], option_name: str
) -> tuple[Choice, ...]:
    """Return string enum option values or raise a clear validation error."""
    return tuple(validate_choice(value, valid, option_name) for value in values)


def existing_path(value: Path | str, option_name: str) -> Path:
    """Return an existing path or raise a clear validation error."""
    path = Path(value)
    if not path.exists():
        raise ValueError(f"{option_name} does not exist: {path}")
    return path
