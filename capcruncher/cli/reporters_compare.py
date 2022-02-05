import itertools
import os
import re
import sys
from typing import Literal, Tuple

import numpy as np
import pandas as pd
from capcruncher.tools.pileup import CoolerBedGraph
from capcruncher.utils import get_cooler_uri
from joblib import Parallel, delayed
from pybedtools import BedTool
from collections import defaultdict


def get_bedgraph_name_from_cooler(cooler_filename):

    filename = os.path.basename(cooler_filename.split(".hdf5")[0])
    viewpoint = cooler_filename.split("::/")[1]
    return f"{filename}_{viewpoint}"


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
    norm_kwargs = {"scale_factor": scale_factor, "regions": normalisation_regions}

    if input_format == "cooler":

        bedgraphs = dict(
            Parallel(n_jobs=n_cores)(
                delayed(
                    lambda fn: (
                        get_bedgraph_name_from_cooler(fn),
                        CoolerBedGraph(fn, region_to_limit=region if region else None)
                        .extract_bedgraph(normalisation=normalisation, **norm_kwargs)
                        .pipe(BedTool.from_dataframe),
                    )
                )(get_cooler_uri(fn, viewpoint, resolution))
                for fn in infiles
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

    return union


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


def get_groups(columns, group_names, group_columns):

    groups = dict()

    for group_name, group_col in zip(group_names, group_columns):
        for col in re.split(r"[,;\s+]", group_col):

            try:
                col = int(col)
                col_name = columns[col]
            except Exception as e:
                col_name = col

            groups[col_name] = group_name

    return groups


# def get_grouped_dataframe(
#     df: pd.DataFrame, group_name: str, groups: dict, agg_func: function
# ):
#     return df.loc[:, groups[group_name]].pipe(agg_func, axis=1)


def summarise(
    infile: os.PathLike,
    output_prefix: os.PathLike = None,
    output_format: Literal["bedgraph", "tsv"] = "bedgraph",
    summary_methods: Tuple[str] = None,
    group_names: Tuple[str] = None,
    group_columns: Tuple[int, str] = None,  # Need to ensure these are 0 based
    suffix: str = "",
    subtraction: bool = False,
):

    df_union = pd.read_csv(infile, sep="\t")
    df_counts = df_union.iloc[:, 3:]
    summary_functions = get_summary_functions(summary_methods)
    groups = (
        get_groups(df_counts.columns, group_names, group_columns)
        if group_names
        else {col: "summary" for col in df_counts.columns}
    )  # Use all columns if no groups provided

    # Perform all groupby aggregations.
    df_agg = (
        df_union.iloc[:, 3:]
        .transpose()  # Transpose to enable groupby funcs
        .groupby(groups)
        .agg([*summary_functions])  # Apply all aggregaions
        .transpose()
        .reset_index()
        .rename(columns={"level_0": "index", "level_1": "aggregation"})
        .set_index("index")
        .join(df_union.iloc[:, :3])
        .fillna(0)
    )

    # Write out groupby aggregations
    for group in group_names:
        if output_format == "bedgraph":
            df_output = df_agg[["chrom", "start", "end", group, "aggregation"]]
            for aggregation, df in df_output.groupby("aggregation"):
                df.drop(columns="aggregation").to_csv(
                    f"{output_prefix}{group}.{aggregation}-summary{suffix}.bedgraph",
                    sep="\t",
                    header=False,
                    index=False,
                )
        elif output_format == "tsv":
            df_output = df_agg[["chrom", "start", "end", "aggregation", group]]
            df_output.to_csv(f"{output_prefix}{group}{suffix}.tsv")

    # Perform permutations
    if subtraction:
        subtractions_performed = list()
        for group_a, group_b in itertools.permutations(group_names, 2):

            subtraction_name = f"{group_a}-{group_b}"
            df_agg[subtraction_name] = df_agg[group_a] - df_agg[group_b]
            subtractions_performed.append(subtraction_name)

        if output_format == "bedgraph":
            for sub in subtractions_performed:
                df_output = df_agg[["chrom", "start", "end", sub, "aggregation"]]
                for aggregation, df in df_output.groupby("aggregation"):
                    df.drop(columns="aggregation").to_csv(
                        f"{output_prefix}{sub}.{aggregation}-subtraction{suffix}.bedgraph",
                        sep="\t",
                        header=False,
                        index=False,
                    )

        elif output_format == "tsv":
            df_output = df_agg[
                [
                    "chrom",
                    "start",
                    "end",
                    "aggregation",
                    *group_names,
                    *subtractions_performed,
                ]
            ]
            df_output.to_csv(f"{output_prefix}subtractions{suffix}.tsv")

    # dataframes_grouped = dict()
    #
    #     for g in [group_a, group_b]:
    #         for func_name, func in summary_functions:
    #             df_grouped = get_grouped_dataframe(df=df_counts,groups=groups, group_name=g, agg_func=func)
    #             df_grouped.to_csv(f"{output_prefix}.{g}.{}")

    # # summary_performed = set()
    # # for group_a, group_b in itertools.permutations(groups, 2):

    # #     df_a = df_counts.loc[:, groups[group_a]]
    # #     df_b = df_counts.loc[:, groups[group_b]]

    # #     counts_for_comparison = list()
    # #     for (group, df) in zip([group_a, group_b], [df_a, df_b]):

    # #         df_summary = pd.concat(
    # #             [
    # #                 pd.Series(df.pipe(summary_functions[summary], axis=1), name=summary)
    # #                 for summary in summary_functions
    # #             ],
    # #             axis=1,
    # #         )

    # #         counts_for_comparison.append(df_summary)

    # #         if not (group, summary_methods) in summary_performed:

    # #             breakpoint()
    # #             if output_format == "bedgraph":

    # #                 for summary in summary_functions:
    # #                     df_bedgraph = pd.concat([df_union.iloc[:, :3], df_summary.loc[:, summary]], axis=1)
    # #                     df_bedgraph.to_csv(f"{output_prefix.replace('.bedgraph', '')}{group}.{summary}.bedgraph", sep="\t", header=False, index=False)

    # #             elif output_format == "tsv":
    # #                 df_bedgraph = pd.concat([df_union.iloc[:, :3], df_summary], axis=1)
    # #                 df_bedgraph.to_csv(f"{output_prefix.replace('.bedgraph', '')}.{group}.summarised.tsv", sep="\t", index=False)

    # #             summary_performed.add((group, summary_methods))

    # #     if subtraction:

    # #         df_subtraction = counts_for_comparison[0] - counts_for_comparison[1]

    # #         if output_format == "bedgraph":

    # #             for summary in summary_functions:
    # #                 df_bedgraph = pd.concat([df_union.iloc[:, :3], df_subtraction.loc[:, summary]], axis=1)
    # #                 df_bedgraph.to_csv(f"{output_prefix.replace('.bedgraph', '').strip('.')}.{group_a}-{group_b}.{summary}-subtraction.bedgraph", sep="\t", header=False, index=False)

    # #         elif output_format == "tsv":
    # #             df_bedgraph = pd.concat([df_union.iloc[:, :3], df_subtraction], axis=1)
    # #             df_bedgraph.to_csv(f"{output_prefix.replace('.bedgraph', '').strip('.')}.{group_a}-{group_b}.{summary}-subtraction.tsv", sep="\t", index=False)
