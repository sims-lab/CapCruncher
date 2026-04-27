from collections import defaultdict
from collections.abc import Callable, Sequence
import itertools
import os
import re
from pathlib import Path

import cooler
from loguru import logger
from pydantic import BaseModel, Field, PositiveFloat, PositiveInt, field_validator, model_validator

import pandas as pd
import polars as pl
from capcruncher.types import (
    CompareFormat,
    Normalisation,
    OutputFormat,
    SummaryMethod,
    VALID_SUMMARY_METHODS,
    existing_path,
    validate_choices,
)

type SummaryFunction = Callable[[pd.Series], float]


class CompareConcatOptions(BaseModel):
    """Validated options for concatenating bedgraphs or cooler pileups."""

    infiles: tuple[Path | str, ...]
    viewpoint: str | None = None
    resolution: int | None = None
    format: CompareFormat = CompareFormat.AUTO
    region: str | None = None
    output: Path | str | None = None
    normalisation: Normalisation = Normalisation.RAW
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
        if self.normalisation == Normalisation.REGION and self.normalisation_regions is None:
            raise ValueError(
                "normalisation_regions is required when normalisation is 'region'."
            )
        if self.normalisation != Normalisation.REGION and self.normalisation_regions is not None:
            raise ValueError(
                "normalisation_regions can only be used when normalisation is 'region'."
            )
        return self


class CompareSummariseOptions(BaseModel):
    """Validated options for summarising concatenated bedgraph tables."""

    infile: Path
    design_matrix: Path | None = None
    output_prefix: Path | str | None = None
    output_format: OutputFormat = OutputFormat.BEDGRAPH
    summary_methods: tuple[SummaryMethod, ...] = (SummaryMethod.MEAN,)
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
            tuple(value or (SummaryMethod.MEAN,)), VALID_SUMMARY_METHODS, "summary_methods"
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


def remove_duplicate_entries(df: pd.DataFrame) -> pd.DataFrame:
    """Removes duplicate coordinates by aggregating values."""
    coordinate_columns = ["chrom", "start", "end"]

    return (
        df.groupby(coordinate_columns)
        .agg("sum")
        .reset_index()
        .sort_values(coordinate_columns)
    )


def concat(
    infiles: Sequence[Path | str],
    viewpoint: str | None = None,
    resolution: int | None = None,
    format: CompareFormat = CompareFormat.AUTO,
    region: str | None = None,
    output: Path | str | None = None,
    normalisation: Normalisation = Normalisation.RAW,
    n_cores: int = 1,
    scale_factor: int = int(1e6),
    normalisation_regions: Path | str | None = None,
) -> dict[str, pd.DataFrame]:
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
    input_format = options.format
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

        def _prepare_bedgraph(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
            df = remove_duplicate_entries(df)
            return df.rename(columns={"score": column_name})

        if input_format == CompareFormat.COOLER:
            from capcruncher.api.pileup import CoolerBedGraph
            from capcruncher.utils import get_cooler_uri

            cooler_uris = [
                get_cooler_uri(fn, viewpoint, options.resolution) for fn in infiles
            ]
            bedgraphs = {
                get_bedgraph_name_from_cooler(uri): _prepare_bedgraph(
                    CoolerBedGraph(
                        uri, region_to_limit=options.region if options.region else None
                    ).extract_bedgraph(
                        normalisation=options.normalisation, **norm_kwargs
                    ),
                    get_bedgraph_name_from_cooler(uri),
                )
                for uri in cooler_uris
            }

        elif input_format == CompareFormat.BEDGRAPH:

            bedgraphs = {
                os.path.basename(fn): _prepare_bedgraph(
                    pd.read_csv(
                        fn,
                        sep="\t",
                        header=None,
                        names=["chrom", "start", "end", "score"],
                    ),
                    os.path.basename(fn),
                )
                for fn in infiles
            }

        else:
            raise NotImplementedError("Auto currently not implemented")

        union = None
        for name, df in bedgraphs.items():
            df = df.rename(columns={"score": name})
            if union is None:
                union = df
            else:
                union = union.merge(df, on=coordinate_columns, how="outer")

        if union is None:
            union = pd.DataFrame(columns=coordinate_columns)

        value_columns = [col for col in union.columns if col not in coordinate_columns]
        if value_columns:
            union[value_columns] = union[value_columns].fillna(0)
        union = union.sort_values(coordinate_columns).reset_index(drop=True)

        if options.output:
            union.to_csv(options.output, sep="\t", index=False)

        union_by_viewpoint[viewpoint] = union

    return union_by_viewpoint


def get_summary_functions(methods: Sequence[str] | None) -> dict[str, SummaryFunction]:
    import numpy as np
    import scipy.stats

    if methods:
        summary_functions = dict()
        for method in methods:
            for package in [np, scipy.stats]:

                if not summary_functions.get(method):
                    try:
                        summary_functions[method] = getattr(package, method)
                    except AttributeError:
                        pass
    else:
        summary_functions = {"mean": getattr(np, "mean")}

    return summary_functions


def get_groups(
    columns: pd.Index | list[str],
    group_names: Sequence[str],
    group_columns: Sequence[str | int],
) -> dict[str, str]:
    """Extracts groups from group_columns and returns a dictionary of column names to group names."""

    groups = dict()

    for group_name, group_col in zip(group_names, group_columns):
        for col in re.split(r"[,;\s+]", group_col):

            try:
                col = int(col)
                col_name = columns[col]
            except Exception:
                col_name = col

            groups[col_name] = group_name

    return groups


def summarise(
    infile: Path | str,
    design_matrix: Path | str | None = None,
    output_prefix: Path | str | None = None,
    output_format: OutputFormat = OutputFormat.BEDGRAPH,
    summary_methods: tuple[SummaryMethod, ...] = (SummaryMethod.MEAN,),
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
        infile=infile,
        design_matrix=design_matrix,
        output_prefix=output_prefix,
        output_format=output_format,
        summary_methods=summary_methods,
        group_names=group_names,
        group_columns=group_columns,
        suffix=suffix,
        perform_subtractions=perform_subtractions,
    )
    logger.info(f"Reading {options.infile}")
    df_union = pd.read_csv(options.infile, sep="\t")
    df_counts = df_union.iloc[:, 3:]

    logger.info("Identifying groups")
    if options.group_columns and options.group_names:
        groups = (
            get_groups(df_counts.columns, options.group_names, options.group_columns)
            if options.group_names
            else {col: "summary" for col in df_counts.columns}
        )  # Use all columns if no groups provided
    
    elif options.design_matrix:
        df_design = pd.read_csv(options.design_matrix, sep=r"\s+|,|\t", engine="python")
        # This design file should look like: sample, condition
        groups = df_design.set_index("sample").to_dict()["condition"]
    else:
        logger.warning("No groups provided, using all columns")
        groups = {col: "summary" for col in df_counts.columns}

    logger.info(f"Extracted groups: {groups}")
    aggregation = defaultdict(list)

    # Invert the groups so conditions are keys
    groups_inverted = defaultdict(list)
    for k, v in groups.items():
        groups_inverted[v].append(k)

    # Convert to polars
    counts = pl.DataFrame(df_counts)
    coordinates = pl.DataFrame(df_union.iloc[:, :3])
    summary_methods = options.summary_methods

    for aggregation_method in summary_methods:
        logger.info(f"Performing aggregation: {aggregation_method}")

        # Apply aggregation method to each group
        for group_name, group in groups_inverted.items():
            colname = f"{group_name}_{aggregation_method}"
            group_counts = getattr(
                counts.select(group), f"{aggregation_method}_horizontal"
            )().alias(colname)
            coordinates = coordinates.with_columns(group_counts)
            aggregation[aggregation_method].append(colname)

        # Perform subtractions
        subtraction = list()
        if perform_subtractions:
            for group_a, group_b in itertools.permutations(groups_inverted, 2):
                group_a_col = f"{group_a}_{aggregation_method}"
                group_b_col = f"{group_b}_{aggregation_method}"

                a = coordinates.select(group_a_col)
                b = coordinates.select(group_b_col)
                diff = a.mean_horizontal() - b.mean_horizontal()
                coordinates = coordinates.with_columns(
                    diff.alias(f"{group_a}-{group_b}")
                )
                subtraction.append(f"{group_a}-{group_b}")

        # Export aggregations
        if options.output_format == OutputFormat.BEDGRAPH:
            # Check that there are no duplicate chrom, start, end coordinates
            coordinates = coordinates.unique(subset=["chrom", "start", "end"])

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
