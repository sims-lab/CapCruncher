"""Identifies differential interactions between conditions."""

import itertools
import os

import pandas as pd
import ray
from loguru import logger
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from capcruncher.api.pileup import cooler_to_bedgraph


@ray.remote
def get_differential_interactions(
    counts: pd.DataFrame,
    design: pd.DataFrame,
    contrast: str,
    group_a: str,
    group_b: str,
    threshold_q: float = 0.05,
    lfc_shrink: bool = False,
):
    """Runs DESeq2 on interaction counts."""
    # Create DeseqDataSet

    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts.T,
        metadata=design,
        design=f"~{contrast}",
        refit_cooks=True,
        inference=inference,
    )

    # Run DESeq2
    dds.deseq2()

    # Get results
    ds = DeseqStats(dds, contrast=[contrast, group_b, group_a], inference=inference)
    df_results = ds.summary()

    if lfc_shrink:
        df_results = ds.lfc_shrink()

    # Filter results
    df_results = df_results.loc[lambda df: df["padj"] <= threshold_q]

    # Sort results
    df_results = (
        df_results.assign(log2FoldChangeAbs=lambda df: df["log2FoldChange"].abs())
        .sort_values(by=["log2FoldChangeAbs", "padj"], ascending=[False, True])
        .drop(columns=["log2FoldChangeAbs"])
    )

    # Add coordinates
    df_results = df_results.assign(
        chrom=lambda df: df.index.str.split(":").str[0],
        start=lambda df: df.index.str.split(":")
        .str[1]
        .str.split("-")
        .str[0]
        .astype(int),
        end=lambda df: df.index.str.split(":").str[1].str.split("-").str[1].astype(int),
    )

    return df_results


def differential(
    interaction_files: list,
    viewpoint: str,
    design_matrix: os.PathLike,
    output_prefix: os.PathLike = "differential_interactions",
    contrast: str = "condition",
    regions_of_interest: os.PathLike = None,
    viewpoint_distance: int = None,
    threshold_count: float = 20,
    threshold_q: float = 0.05,
):
    """Identifies differential interactions between conditions.

    Parses a list of cooler files containg reporter counts from at least two conditions with
    two or more replicates for a single capture probe and outputs differential interaction
    results. Following filtering to ensure that the number of interactions is above the required
    threshold, PyDeseq2 is used to run a compatison after
    fitting a negative binomial model to the interaction counts.The options to filter
    results can be filtered by a minimum mean value (threshold_mean) and/or
    maximum q-value (threshold-q) are also provided.


    Args:
        interaction_files (list): List of cooler files.
        viewpoint (str): Name of capture probe. MUST match one viewpoint within the HDF5 files.
        design_matrix (os.PathLike): Design matrix to use for grouping samples. (N_SAMPLES * METADATA).
        output_prefix (os.PathLike, optional): Output prefix for differntial interactions. Defaults to 'differential'.
        contrast (str, optional): Column to use for grouping. Defaults to 'condition'.
        regions_of_interest (os.PathLike, optional): BED file of regions of interest. Defaults to None.
        viewpoint_distance (int, optional): Distance from viewpoint to include. Defaults to 500_000.
        threshold_count (float, optional): Minimum number of reported interactions required. Defaults to 20.
        threshold_q (float, optional): Maximum q-value for output. Defaults to 0.05.
        threshold_mean (float, optional): Minimum mean value for output. Defaults to 0.
    """
    # Load design matrix
    logger.info("Loading design matrix.")
    df_design = pd.read_table(
        design_matrix, index_col=0, sep=r"\s+|,|\t", engine="python"
    )

    # Set-up tasks for bedgraph generation
    logger.info("Validating viewpoint distance and regions of interest.")
    assert len(interaction_files) >= 2, "No interaction files provided."
    assert regions_of_interest or viewpoint_distance, (
        "No regions of interest or viewpoint distance provided."
    )

    logger.info("Extracting interaction counts.")

    if regions_of_interest:
        logger.info(
            f"Using supplied regions of interest file {regions_of_interest} to restrict analysis"
        )
    else:
        logger.info(
            f"Using distance from viewpoint of {viewpoint_distance} to restrict analysis"
        )

    bedgraph_futures = dict()
    for interaction_file in interaction_files:
        file_name = os.path.basename(interaction_file.replace(".hdf5", ""))
        future = cooler_to_bedgraph.remote(
            clr=f"{interaction_file}::{viewpoint}",
            regions_of_interest=regions_of_interest,
            viewpoint_distance=viewpoint_distance,
        )
        bedgraph_futures[file_name] = future

    # Execute tasks
    bedgraphs = {k: ray.get(v) for k, v in bedgraph_futures.items()}

    logger.info("Concatenating interactions.")
    # Concatenate bedgraphs
    df_counts = pd.concat(
        [
            bg.assign(
                coord=lambda df: df["chrom"].astype(str)
                + ":"
                + df["start"].astype(str)
                + "-"
                + df["end"].astype(str)
            )
            .set_index("coord")
            .drop(columns=["chrom", "start", "end"])
            .rename(columns={"count": name})
            for name, bg in bedgraphs.items()
        ],
        axis=1,
    ).fillna(0)

    # Filter out any interacting fragments with less than threshold_counts
    logger.info(f"Removing interactions with less than {threshold_count} counts.")
    df_counts = df_counts.loc[lambda df: (df >= threshold_count).all(axis=1)]

    # At the time of writing. PyDeseq2 doese not support multiple comparisons.
    # Therefore, we need to run a separate DESeq2 analysis for each comparison.

    # Get all comparisons
    possible_contrasts = df_design[contrast].unique()
    comparisons = list(itertools.combinations(possible_contrasts, 2))

    # Run comparisons
    comparison_futures = dict()
    for group_a, group_b in comparisons:
        # Filter design matrix
        df_design_sub = df_design.loc[lambda df: df[contrast].isin([group_a, group_b])]

        # Filter counts
        df_counts_sub = df_counts.loc[:, df_design_sub.index]

        # Get differential interactions
        result = get_differential_interactions.remote(
            df_counts_sub,
            df_design_sub,
            contrast,
            threshold_q=threshold_q,
            group_a=group_a,
            group_b=group_b,
        )

        comparison_futures[(group_a, group_b)] = result

    # Execute tasks
    for (group_a, group_b), future in comparison_futures.items():
        logger.info(f"Running comparison: {group_a} vs {group_b}")
        df_results = ray.get(future)

        # Write result
        df_results.to_csv(
            f"{output_prefix}.{group_a}_vs_{group_b}.csv",
            sep=",",
            index=True,
            header=True,
        )
