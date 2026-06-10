import itertools
import os
import re
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import cast

import cooler
import polars as pl
from loguru import logger
from pydantic import (
    BaseModel,
    PositiveFloat,
    PositiveInt,
    field_validator,
    model_validator,
)

from capcruncher.types import (
    VALID_SUMMARY_METHODS,
    CompareFormat,
    Normalisation,
    OutputFormat,
    SummaryMethod,
    existing_path,
    validate_choices,
)


class CompareConcatOptions(BaseModel):
    """Validated options for concatenating bedgraphs or cooler pileups."""

    infiles: tuple[Path | str, ...]
    viewpoint: str | None = None
    resolution: int | None = None
    format: CompareFormat | str = CompareFormat.AUTO
    region: str | None = None
    output: Path | str | None = None
    normalisation: Normalisation | str = Normalisation.RAW
    n_cores: PositiveInt = 1
    scale_factor: PositiveFloat = 1e6
    normalisation_regions: Path | str | None = None

    @field_validator("infiles", mode="before")
    @classmethod
    def validate_infiles(cls, value: Sequence[Path | str]) -> tuple[Path | str, ...]:
        values = tuple(value)
        if not values:
            raise ValueError("At least one input file is required.")
        return values

    @field_validator("normalisation_regions", mode="before")
    @classmethod
    def empty_region_to_none(cls, value: Path | str | None) -> Path | str | None:
        return None if value == "" else value

    @model_validator(mode="after")
    def validate_normalisation_regions(self) -> "CompareConcatOptions":
        if (
            self.normalisation == Normalisation.REGION
            and self.normalisation_regions is None
        ):
            raise ValueError(
                "normalisation_regions is required when normalisation is 'region'."
            )
        if (
            self.normalisation != Normalisation.REGION
            and self.normalisation_regions is not None
        ):
            raise ValueError(
                "normalisation_regions can only be used when normalisation is 'region'."
            )
        return self


class CompareSummariseOptions(BaseModel):
    """Validated options for summarising concatenated bedgraph tables."""

    infile: Path
    design_matrix: Path | None = None
    output_prefix: Path | str | None = None
    output_format: OutputFormat | str = OutputFormat.BEDGRAPH
    summary_methods: tuple[SummaryMethod | str, ...] = (SummaryMethod.MEAN,)
    group_names: tuple[str, ...] = ()
    group_columns: tuple[str | int, ...] = ()
    suffix: str = ""
    perform_subtractions: bool = False

    @field_validator("infile", mode="before")
    @classmethod
    def validate_infile(cls, value: Path | str) -> Path:
        return existing_path(value, "infile")

    @field_validator("design_matrix", mode="before")
    @classmethod
    def validate_design_matrix(cls, value: Path | str | None) -> Path | None:
        if value in (None, ""):
            return None
        return existing_path(value, "design_matrix")

    @field_validator("summary_methods", mode="before")
    @classmethod
    def validate_summary_methods(
        cls, value: Sequence[str] | None
    ) -> tuple[SummaryMethod, ...]:
        return validate_choices(
            tuple(value or (SummaryMethod.MEAN,)),
            VALID_SUMMARY_METHODS,
            "summary_methods",
        )

    @field_validator("group_names", mode="before")
    @classmethod
    def validate_group_names(cls, value: Sequence[str] | None) -> tuple[str, ...]:
        return tuple(value or ())

    @field_validator("group_columns", mode="before")
    @classmethod
    def validate_group_columns(
        cls, value: Sequence[str | int] | None
    ) -> tuple[str | int, ...]:
        return tuple(value or ())


def get_bedgraph_name_from_cooler(cooler_filename: str) -> str:
    filename = os.path.basename(cooler_filename.split(".hdf5")[0])
    viewpoint = cooler_filename.split("::/")[1]
    return f"{filename}_{viewpoint}"


def _to_polars(df) -> pl.DataFrame:
    return df if isinstance(df, pl.DataFrame) else pl.from_pandas(df)


def remove_duplicate_entries(df) -> pl.DataFrame:
    """Removes duplicate coordinates by aggregating values."""
    coordinate_columns = ["chrom", "start", "end"]

    return (
        _to_polars(df)
        .group_by(coordinate_columns)
        .agg(pl.exclude(coordinate_columns).sum())
        .sort(coordinate_columns)
    )


def concat(
    infiles: Sequence[Path | str],
    viewpoint: str | None = None,
    resolution: int | None = None,
    format: CompareFormat | str = CompareFormat.AUTO,
    region: str | None = None,
    output: Path | str | None = None,
    normalisation: Normalisation | str = Normalisation.RAW,
    n_cores: int = 1,
    scale_factor: int = int(1e6),
    normalisation_regions: Path | str | None = None,
) -> dict[str, pl.DataFrame]:
    """Concatenate bedgraphs or cooler-derived bedgraphs by viewpoint."""
    options = CompareConcatOptions(
        infiles=tuple(infiles),
        viewpoint=viewpoint,
        resolution=resolution,
        format=format,
        region=region,
        output=output,
        normalisation=normalisation,
        n_cores=n_cores,
        scale_factor=scale_factor,
        normalisation_regions=normalisation_regions,
    )
    input_format = cast(CompareFormat, options.format)
    normalisation = cast(Normalisation, options.normalisation)
    norm_kwargs = {
        "scale_factor": options.scale_factor,
        "region": options.normalisation_regions,
    }
    infiles = [os.fspath(infile) for infile in options.infiles]

    if not options.viewpoint:
        viewpoints = [vp.strip("/") for vp in cooler.fileops.list_coolers(infiles[0])]
    else:
        viewpoints = [
            options.viewpoint,
        ]

    union_by_viewpoint = dict()

    for viewpoint in viewpoints:
        coordinate_columns = ["chrom", "start", "end"]

        def _prepare_bedgraph(df, column_name: str) -> pl.DataFrame:
            df = remove_duplicate_entries(df)
            rename_map = {
                column: column_name
                for column in ("score", "count")
                if column in df.columns
            }
            return df.rename(rename_map)

        if input_format == CompareFormat.COOLER:
            from capcruncher.api.interactions.bedgraph import CoolerBedGraph
            from capcruncher.utils import get_cooler_uri

            cooler_uris = [
                get_cooler_uri(fn, viewpoint, options.resolution) for fn in infiles
            ]
            bedgraphs = {
                get_bedgraph_name_from_cooler(uri): _prepare_bedgraph(
                    CoolerBedGraph(
                        uri, region_to_limit=options.region if options.region else None
                    ).extract_bedgraph(normalisation=normalisation, **norm_kwargs),
                    get_bedgraph_name_from_cooler(uri),
                )
                for uri in cooler_uris
            }

        elif input_format == CompareFormat.BEDGRAPH:
            bedgraphs = {
                os.path.basename(fn): _prepare_bedgraph(
                    pl.read_csv(
                        fn,
                        separator="\t",
                        has_header=False,
                        new_columns=["chrom", "start", "end", "score"],
                    ),
                    os.path.basename(fn),
                )
                for fn in infiles
            }

        else:
            raise NotImplementedError("Auto currently not implemented")

        union = None
        for _, df in bedgraphs.items():
            if union is None:
                union = df
            else:
                union = union.join(df, on=coordinate_columns, how="full", coalesce=True)

        if union is None:
            union = pl.DataFrame(
                schema={column: pl.String for column in coordinate_columns}
            )

        value_columns = [col for col in union.columns if col not in coordinate_columns]
        if value_columns:
            union = union.with_columns(pl.col(value_columns).fill_null(0))
        union = union.sort(coordinate_columns)

        if options.output:
            union.write_csv(options.output, separator="\t")

        union_by_viewpoint[viewpoint] = union

    return union_by_viewpoint


def get_groups(
    columns: list[str],
    group_names: Sequence[str],
    group_columns: Sequence[str | int],
) -> dict[str, str]:
    """Extracts groups from group_columns and returns a dictionary of column names to group names."""

    groups = dict()

    for group_name, group_col in zip(group_names, group_columns, strict=True):
        for col in re.split(r"[,;\s+]", str(group_col)):
            try:
                col = int(col)
                col_name = columns[col]
            except Exception:
                col_name = col

            groups[col_name] = group_name

    return groups


def _read_design_matrix(path: Path) -> pl.DataFrame:
    lines = path.read_text().splitlines()
    rows = [re.split(r"\s+|,|\t", line.strip()) for line in lines if line.strip()]
    if not rows:
        return pl.DataFrame()
    header, *body = rows
    return pl.DataFrame([dict(zip(header, row, strict=False)) for row in body])


def summarise(
    infile: Path | str,
    design_matrix: Path | str | None = None,
    output_prefix: Path | str | None = None,
    output_format: OutputFormat | str = OutputFormat.BEDGRAPH,
    summary_methods: Sequence[SummaryMethod | str] = (SummaryMethod.MEAN,),
    group_names: tuple[str, ...] | None = None,
    group_columns: tuple[str | int, ...] | None = None,
    suffix: str = "",
    perform_subtractions: bool = False,
) -> None:
    """Summarise a concatenated bedgraph table by group.

    Only ``mean`` summaries are currently supported. Unsupported methods raise
    ``ValueError`` before data processing.
    """
    options = CompareSummariseOptions(
        infile=Path(infile),
        design_matrix=Path(design_matrix) if design_matrix is not None else None,
        output_prefix=output_prefix,
        output_format=output_format,
        summary_methods=tuple(summary_methods),
        group_names=tuple(group_names or ()),
        group_columns=tuple(group_columns or ()),
        suffix=suffix,
        perform_subtractions=perform_subtractions,
    )
    logger.info(f"Reading {options.infile}")
    _header = pl.read_csv(options.infile, separator="\t", n_rows=0).columns
    df_union = pl.read_csv(
        options.infile,
        separator="\t",
        schema_overrides={col: pl.Float64 for col in _header[3:]},
    )
    count_columns = df_union.columns[3:]

    logger.info("Identifying groups")
    if options.group_columns and options.group_names:
        groups = (
            get_groups(count_columns, options.group_names, options.group_columns)
            if options.group_names
            else {col: "summary" for col in count_columns}
        )  # Use all columns if no groups provided

    elif options.design_matrix:
        df_design = _read_design_matrix(options.design_matrix)
        # This design file should look like: sample, condition
        groups = dict(
            zip(
                df_design.get_column("sample").to_list(),
                df_design.get_column("condition").to_list(),
                strict=False,
            )
        )
    else:
        logger.warning("No groups provided, using all columns")
        groups = {col: "summary" for col in count_columns}

    logger.info(f"Extracted groups: {groups}")
    aggregation = defaultdict(list)

    # Invert the groups so conditions are keys
    groups_inverted = defaultdict(list)
    for k, v in groups.items():
        groups_inverted[v].append(k)

    counts = df_union.select(count_columns)
    coordinates = df_union.select(df_union.columns[:3])
    summary_methods = options.summary_methods

    for aggregation_method in summary_methods:
        logger.info(f"Performing aggregation: {aggregation_method}")

        # Apply aggregation method to each group
        for group_name, group in groups_inverted.items():
            colname = f"{group_name}_{aggregation_method}"
            group_counts = counts.select(
                pl.mean_horizontal(group).alias(colname)
            ).get_column(colname)
            coordinates = coordinates.with_columns(group_counts)
            aggregation[aggregation_method].append(colname)

        # Perform subtractions
        subtraction = list()
        if perform_subtractions:
            for group_a, group_b in itertools.permutations(groups_inverted, 2):
                group_a_col = f"{group_a}_{aggregation_method}"
                group_b_col = f"{group_b}_{aggregation_method}"

                coordinates = coordinates.with_columns(
                    (pl.col(group_a_col) - pl.col(group_b_col)).alias(
                        f"{group_a}_vs_{group_b}"
                    )
                )
                subtraction.append(f"{group_a}_vs_{group_b}")

        # Export aggregations
        if options.output_format == OutputFormat.BEDGRAPH:
            # Check that there are no duplicate chrom, start, end coordinates
            coordinates = coordinates.unique(
                subset=["chrom", "start", "end"], maintain_order=True
            ).sort(["chrom", "start", "end"])

            # Write the output
            for aggregation_method, group_names in aggregation.items():
                for group_name in group_names:
                    df_output = coordinates.select(
                        ["chrom", "start", "end", group_name]
                    )

                    group_name_cleaned = re.sub(
                        "|".join([*summary_methods, "_"]), "", group_name
                    )
                    outfile = f"{options.output_prefix}{group_name_cleaned}.{aggregation_method}-summary{options.suffix}.bedgraph"

                    logger.info(
                        f"Writing {group_name} {aggregation_method} to {outfile}"
                    )
                    df_output.write_csv(outfile, separator="\t", include_header=False)

            for sub in subtraction:
                df_output = coordinates.select(["chrom", "start", "end", sub])
                outfile = f"{options.output_prefix}{sub}.{aggregation_method}-subtraction{options.suffix}.bedgraph"
                logger.info(f"Writing {sub} {aggregation_method} to {outfile}")
                df_output.write_csv(outfile, separator="\t", include_header=False)
        elif options.output_format == OutputFormat.TSV:
            df_output = coordinates
            df_output.write_csv(
                f"{options.output_prefix}{options.suffix}.tsv",
                separator="\t",
                include_header=True,
            )
