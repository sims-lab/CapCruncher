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

    groups = defaultdict(list)

    for name, col in zip(group_names, group_columns):
        for c in re.split(r"[,;\s+]", col):

            try:
                c = int(c)
                col_name = columns[c]
            except Exception as e:
                col_name = c

            groups[name].append(col_name)

    return groups


def summarise(
    infile: os.PathLike,
    output_prefix: os.PathLike = None,
    output_format: Literal["bedgraph", "tsv"] = "bedgraph",
    summary_methods: Tuple[str] = None,
    group_names: Tuple[str] = None,
    group_columns: Tuple[int] = None,  # Need to ensure these are 0 based
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

    summary_performed = set()
    for group_a, group_b in itertools.permutations(groups, 2):

        df_a = df_counts.loc[:, groups[group_a]]
        df_b = df_counts.loc[:, groups[group_b]]

        counts_for_comparison = list()
        for (group, df) in zip([group_a, group_b], [df_a, df_b]):

            df_summary = pd.concat(
                [
                    pd.Series(df.pipe(summary_functions[summary], axis=1), name=summary)
                    for summary in summary_functions
                ],
                axis=1,
            )

            counts_for_comparison.append(df_summary)

            if not (group, summary_methods) in summary_performed:
                if output_format == "bedgraph":

                    for summary in summary_functions:
                        df_bedgraph = pd.concat([df_union.iloc[:, :3], df_summary.loc[:, summary]], axis=1)
                        df_bedgraph.to_csv(f"{output_prefix.replace('.bedgraph', '')}{group}.{summary}.bedgraph", sep="\t", header=False, index=False)
                
                elif output_format == "tsv":
                    df_bedgraph = pd.concat([df_union.iloc[:, :3], df_summary], axis=1)
                    df_bedgraph.to_csv(f"{output_prefix.replace('.bedgraph', '')}.{group}.summarised.tsv", sep="\t", index=False)
                
                summary_performed.add((group, summary_methods))


        if subtraction:

            df_subtraction = counts_for_comparison[0] - counts_for_comparison[1]
            
            if output_format == "bedgraph":

                for summary in summary_functions:
                    df_bedgraph = pd.concat([df_union.iloc[:, :3], df_subtraction.loc[:, summary]], axis=1)
                    df_bedgraph.to_csv(f"{output_prefix.replace('.bedgraph', '').strip('.')}.{group_a}-{group_b}.{summary}-subtraction.bedgraph", sep="\t", header=False, index=False)
                
            elif output_format == "tsv":
                df_bedgraph = pd.concat([df_union.iloc[:, :3], df_subtraction], axis=1)
                df_bedgraph.to_csv(f"{output_prefix.replace('.bedgraph', '').strip('.')}.{group_a}-{group_b}.{summary}-subtraction.tsv", sep="\t", index=False)
