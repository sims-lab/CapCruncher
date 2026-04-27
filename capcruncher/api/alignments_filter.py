import pathlib
import tempfile
from pathlib import Path
from typing import cast

import pandas as pd
import polars as pl
from loguru import logger
from pydantic import BaseModel, field_validator

from capcruncher.api.filter import CCSliceFilter, TiledCSliceFilter, TriCSliceFilter
from capcruncher.api.io import parse_bam
from capcruncher.api.statistics import SliceFilterStatsList
from capcruncher.types import (
    Assay,
    ReadType,
    VALID_ASSAYS,
    VALID_READ_TYPES,
    existing_path,
    validate_choice,
)

SLICE_FILTERS = {
    "capture": CCSliceFilter,
    "tri": TriCSliceFilter,
    "tiled": TiledCSliceFilter,
}


class AlignmentFilterOptions(BaseModel):
    """Validated options for alignment filtering."""

    bam: Path | str
    annotations: Path | str
    custom_filtering: Path | str | None = None
    output_prefix: Path | str = "reporters"
    statistics: Path = Path("filtering_stats.json")
    method: Assay | str = Assay.CAPTURE
    sample_name: str | None = ""
    read_type: ReadType | str = ReadType.FLASHED
    fragments: bool = True

    @field_validator("bam", "annotations", mode="before")
    @classmethod
    def validate_required_paths(cls, value: Path | str, info) -> Path:
        return existing_path(value, info.field_name)

    @field_validator("custom_filtering", mode="before")
    @classmethod
    def validate_custom_filtering(cls, value: Path | str | None) -> Path | None:
        if value in (None, ""):
            return None
        return existing_path(value, "custom_filtering")

    @field_validator("method", mode="before")
    @classmethod
    def validate_method(cls, value: str | Assay) -> Assay:
        return validate_choice(value, VALID_ASSAYS, "method")

    @field_validator("read_type", mode="before")
    @classmethod
    def validate_read_type(cls, value: str | ReadType) -> ReadType:
        return validate_choice(value, VALID_READ_TYPES, "read_type")


def remove_unused_categories(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for column in df.select_dtypes(include="category").columns:
        df[column] = df[column].cat.remove_unused_categories()
    return df


def merge_annotations(slices: Path | str, annotations: Path | str) -> pd.DataFrame:
    """
    Merges a parquet file containing slice information with a parquet file containing
    annotation information.

    Args:
        slices (os.PathLike): Path to parquet file containing slice information
        annotations (os.PathLike): Path to parquet file containing annotation information

    Returns:
        pd.DataFrame: Merged dataframe
    """

    logger.info("Opening annotations")

    with pl.StringCache():
        join_key_types = {
            "slice_name": pl.Utf8,
            "chrom": pl.Utf8,
            "start": pl.Int64,
        }

        df_slices = pl.scan_parquet(slices).with_columns(
            pl.col(column).cast(dtype) for column, dtype in join_key_types.items()
        )
        df_annotations = pl.scan_parquet(annotations).rename(
            {"Chromosome": "chrom", "Start": "start", "End": "end"}
        ).with_columns(
            pl.col(column).cast(dtype) for column, dtype in join_key_types.items()
        )

        df_slices = df_slices.join(
            df_annotations, on=["slice_name", "chrom", "start"], how="inner"
        )
        df_slices = df_slices.unique(subset=["slice_name"])

        return df_slices.collect().to_pandas()


def filter(
    bam: Path | str,
    annotations: Path | str,
    custom_filtering: Path | str | None = None,
    output_prefix: Path | str = "reporters",
    statistics: Path | str = "",
    method: Assay | str = Assay.CAPTURE,
    sample_name: str | None = "",
    read_type: ReadType | str = ReadType.FLASHED,
    fragments: bool = True,
) -> None:
    """Remove unwanted aligned slices and identify reporters.

    Parses a BAM file, joins it to annotation parquet output, and applies the
    filter set for ``capture``, ``tri``, or ``tiled`` assays.

    Common filters include:

    - Removal of unmapped slices
    - Removal of excluded/blacklisted slices
    - Removal of non-capture fragments
    - Removal of multi-capture fragments
    - Removal of non-reporter fragments
    - Removal of fragments with duplicated coordinates.

    For specific filtering for each of the three methods see:

    - :class:`CCSliceFilter <capcruncher.tools.filter.CCSliceFilter>`
    - :class:`TriCSliceFilter <capcruncher.tools.filter.TriCSliceFilter>`
    - :class:`TiledCSliceFilter <capcruncher.tools.filter.TiledCSliceFilter>`


    In addition to outputting valid reporter fragments and slices separated by capture probe,
    this script also provides statistics on the number of read/slices filtered at each stage,
    and the number of cis/trans reporters for each probe.

    Notes:

     Whilst the script is capable of processing any annotations in tsv format, provided
     that the correct columns are present. It is highly recomended that the "annotate"
     subcomand is used to generate this file.

     Slice filtering is currently hard coded into each filtering class. This will be
     modified in a future update to enable custom filtering orders.


    Args:
        bam: Input BAM file.
        annotations: Annotation parquet generated by ``alignments annotate``.
        custom_filtering: Optional YAML file defining filter stage order.
        output_prefix: Prefix for reporter parquet outputs.
        statistics: Output path for JSON filter statistics.
        method: Assay filter to use: ``capture``, ``tri``, or ``tiled``.
        sample_name: Sample name written to statistics.
        read_type: Read type written to statistics: ``flashed`` or ``pe``.
        fragments: Whether to write fragment-level reporter parquet.

    Raises:
        ValueError: If user-facing option values or required paths are invalid.
    """
    options = AlignmentFilterOptions(
        bam=bam,
        annotations=annotations,
        custom_filtering=custom_filtering,
        output_prefix=output_prefix,
        statistics=Path(statistics) if statistics else Path("filtering_stats.json"),
        method=method,
        sample_name=sample_name or "",
        read_type=read_type,
        fragments=fragments,
    )

    with logger.catch():
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir) / "tmp.parquet"

            logger.info("Loading bam file")
            # Its faster to write to parquet and then read it back than to join both dataframes with pandas
            parse_bam(options.bam).to_parquet(tmp)

            # Join bam file with annotations
            logger.info("Merging bam file with annotations")
            df_alignment = merge_annotations(tmp, options.annotations)

            # Make sure that the blacklist column is present
            if "blacklist" not in df_alignment.columns:
                df_alignment["blacklist"] = 0

        # Initialise SliceFilter
        # If no custom filtering, will use the class default.
        method = cast(Assay, options.method)
        read_type = cast(ReadType, options.read_type)
        slice_filter_class = SLICE_FILTERS[method.value]
        slice_filter = slice_filter_class(
            slices=df_alignment,
            sample_name=options.sample_name,
            read_type=read_type.value,
            filter_stages=options.custom_filtering,
        )

        # Filter slices using the slice_filter
        logger.info(f"Filtering slices with method: {method}")
        slice_filter.filter_slices()

        # Extract statistics
        logger.info("Extracting statistics")
        stats_list = SliceFilterStatsList.from_list(slice_filter.filtering_stats)
        with open(options.statistics, "w") as f:
            f.write(stats_list.model_dump_json())

        # Write output
        df_slices = slice_filter.slices
        df_slices_with_viewpoint = slice_filter.slices_with_viewpoint
        df_capture = slice_filter.captures

        if fragments:
            logger.info("Writing reporters at the fragment level")
            df_fragments = (
                slice_filter_class(df_slices)
                .fragments.join(
                    df_capture["capture"], lsuffix="_slices", rsuffix="_capture"
                )
                .rename(
                    columns={
                        "capture_slices": "capture",
                        "capture_capture": "viewpoint",
                    }
                )
            )
            df_fragments = remove_unused_categories(df_fragments)

            df_fragments.to_parquet(
                f"{options.output_prefix}.fragments.parquet",
                compression="snappy",
                engine="pyarrow",
            )

        logger.info("Writing reporters slices")

        # Enforce dtype for parent_id
        df_slices_with_viewpoint = df_slices_with_viewpoint.assign(
            parent_id=lambda df: df["parent_id"].astype("int64")
        ).drop_duplicates("slice_id")

        # Convert objects to category
        to_convert = df_slices_with_viewpoint.select_dtypes(include="object").columns
        df_slices_with_viewpoint[to_convert] = df_slices_with_viewpoint[
            to_convert
        ].astype("category")
        df_slices_with_viewpoint = remove_unused_categories(df_slices_with_viewpoint)

        df_slices_with_viewpoint.to_parquet(
            f"{options.output_prefix}.slices.parquet",
            compression="snappy",
            engine="pyarrow",
        )

        logger.info("Completed analysis of BAM file")
