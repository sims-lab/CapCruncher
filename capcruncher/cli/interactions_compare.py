import itertools
from loguru import logger
import os
import re
from typing import Literal, Tuple, List, Union, Dict
import cooler

import pandas as pd
import polars as pl


from capcruncher.api.pileup import CoolerBedGraph
from capcruncher.utils import get_cooler_uri
from joblib import Parallel, delayed
from pybedtools import BedTool
from collections import defaultdict


def get_bedgraph_name_from_cooler(cooler_filename):

    filename = os.path.basename(cooler_filename.split(".hdf5")[0])
    viewpoint = cooler_filename.split("::/")[1]
    return f"{filename}_{viewpoint}"


def remove_duplicate_entries(df: pd.DataFrame) -> pd.DataFrame:
    """Removes duplicate coordinates by aggregating values."""

    return (
        df.groupby(["chrom", "start", "end"])
        .agg("sum")
        .reset_index()
        .sort_values(["chrom", "start", "end"])
    )


def concat(
    infiles: Tuple[os.PathLike],
    viewpoint: str = None,
    resolution: int = None,
    format: Literal["auto", "cooler", "bedgraph"] = "auto",
    region: str = None,
    output: os.PathLike = None,
    normalisation: Literal["raw", "n_cis", "region"] = "raw",
    n_cores: int = 1,
    scale_factor: int = int(1e6),
    normalisation_regions: os.PathLike = None,
):

    input_format = format
    norm_kwargs = {"scale_factor": scale_factor, "region": normalisation_regions}

    if not viewpoint:
        viewpoints = [vp.strip("/") for vp in cooler.fileops.list_coolers(infiles[0])]
    else:
        viewpoints = [
            viewpoint,
        ]

    union_by_viewpoint = dict()

    for viewpoint in viewpoints:

        if input_format == "cooler":

            cooler_uris = [get_cooler_uri(fn, viewpoint, resolution) for fn in infiles]
            bedgraphs = dict(
                Parallel(n_jobs=n_cores)(
                    delayed(
                        lambda uri: (
                            get_bedgraph_name_from_cooler(uri),
                            CoolerBedGraph(
                                uri, region_to_limit=region if region else None
                            )
                            .extract_bedgraph(
                                normalisation=normalisation, **norm_kwargs
                            )
                            .pipe(BedTool.from_dataframe),
                        )
                    )(uri)
                    for uri in cooler_uris
                )
            )

        elif input_format == "bedgraph":

            bedgraphs = {os.path.basename(fn): BedTool(fn) for fn in infiles}

        else:
            raise NotImplementedError("Auto currently not implemented")

        union = (
            BedTool()
            .union_bedgraphs(i=[bt.fn for bt in bedgraphs.values()])
            .to_dataframe(
                disable_auto_names=True,
                names=["chrom", "start", "end", *list(bedgraphs.keys())],
            )
        )

        if output:
            union.to_csv(output, sep="\t", index=False)

        union_by_viewpoint[viewpoint] = union

    return union_by_viewpoint


def get_summary_functions(methods):
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
    columns: Union[pd.Index, list],
    group_names: List[str],
    group_columns: List[Union[str, int]],
) -> Dict[str, str]:
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
    infile: os.PathLike,
    design_matrix: os.PathLike = None,
    output_prefix: os.PathLike = None,
    output_format: Literal["bedgraph", "tsv"] = "bedgraph",
    summary_methods: Tuple[Literal['mean']] = ("mean",),
    group_names: Tuple[str] = None,
    group_columns: Tuple[int, str] = None,  # Need to ensure these are 0 based
    suffix: str = "",
    perform_subtractions: bool = False,
):

    logger.info(f"Reading {infile}")
    df_union = pd.read_csv(infile, sep="\t")
    df_counts = df_union.iloc[:, 3:]

    logger.info("Identifying groups")
    if group_columns and group_names:
        groups = (
            get_groups(df_counts.columns, group_names, group_columns)
            if group_names
            else {col: "summary" for col in df_counts.columns}
        )  # Use all columns if no groups provided
    
    elif design_matrix:
        df_design = pd.read_csv(design_matrix, sep=r"\s+|,|\t", engine="python")
        # This design file should look like: sample, condition
        groups = df_design.set_index("sample").to_dict()["condition"]
    else:
        logger.warning("No groups provided, using all columns")

    logger.info(f"Extracted groups: {groups}")
    aggregation = defaultdict(list)
    subtraction = list()

    # Invert the groups so conditions are keys
    groups_inverted = defaultdict(list) 
    for k, v in groups.items():
        groups_inverted[v].append(k)

    # Convert to polars
    counts = pl.DataFrame(df_counts)
    coordinates = pl.DataFrame(df_union.iloc[:, :3])
    summary_methods = ['mean', ] if not summary_methods else summary_methods

    for aggregation_method in summary_methods:
        
        assert aggregation_method in ["mean"], f"Invalid aggregation method {aggregation_method}"
        logger.info(f"Performing aggregation: {aggregation_method}")


        # Apply aggregation method to each group
        for group_name, group in groups_inverted.items():

            colname = f'{group_name}_{aggregation_method}'
            group_counts = getattr(counts.select(group), f'{aggregation_method}_horizontal')().alias(colname)
            coordinates = coordinates.with_columns(group_counts)
            aggregation[aggregation_method].append(colname)
        
        # Perform subtractions
        subtraction = list()
        if perform_subtractions:
            for group_a, group_b in itertools.permutations(groups_inverted, 2):

                group_a_col = f'{group_a}_{aggregation_method}'
                group_b_col = f'{group_b}_{aggregation_method}'

                a = coordinates.select(group_a_col)
                b = coordinates.select(group_b_col)
                diff = a.mean_horizontal() - b.mean_horizontal()
                coordinates = coordinates.with_columns(diff.alias(f"{group_a}-{group_b}"))
                subtraction.append(f"{group_a}-{group_b}")

        # Export aggregations
        if output_format == "bedgraph":
            
            # Check that there are no duplicate chrom, start, end coordinates
            coordinates =  coordinates.unique(subset=["chrom", "start", "end"])

            # Write the output
            for aggregation_method, group_names in aggregation.items():
                for group_name in group_names:
                    df_output = coordinates.select(["chrom", "start", "end", group_name])
                    
                    group_name_cleaned = re.sub('|'.join([*summary_methods, '_']), '', group_name) # Remove the aggregation method from the group name
                    outfile = f"{output_prefix}{group_name_cleaned}.{aggregation_method}-summary{suffix}.bedgraph"
                    
                    logger.info(f"Writing {group_name} {aggregation_method} to {outfile}")
                    df_output.write_csv(outfile, separator="\t", include_header=False)
            
            for sub in subtraction:
                df_output = coordinates.select(["chrom", "start", "end", sub])
                outfile = f"{output_prefix}{sub}.{aggregation_method}-subtraction{suffix}.bedgraph"
                logger.info(f"Writing {sub} {aggregation_method} to {outfile}")
                df_output.write_csv(outfile, separator="\t", include_header=False)
        
        elif output_format == "tsv":
            df_output = coordinates
            df_output.write_csv(f"{output_prefix}{suffix}.tsv", separator="\t", include_header=True)
        