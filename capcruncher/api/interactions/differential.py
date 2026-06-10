"""Identifies differential interactions between conditions."""

import itertools
import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from capcruncher.api.interactions.bedgraph import cooler_to_bedgraph


def _load_pydeseq2() -> tuple[Any, Any, Any]:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
        from pydeseq2.ds import DeseqStats
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.startswith("pydeseq2"):
            raise ModuleNotFoundError(
                "PyDESeq2 is required for differential interactions. "
                "Install CapCruncher with the 'differential' extra."
            ) from exc
        raise

    return DeseqDataSet, DefaultInference, DeseqStats


def _results_dataframe(deseq_stats: Any) -> pd.DataFrame:
    """Return the current PyDESeq2 results table.

    PyDESeq2 0.5 stores results on ``results_df`` and its ``summary`` and
    ``lfc_shrink`` methods mutate that table in place.
    """
    results = getattr(deseq_stats, "results_df", None)
    if results is None:
        raise RuntimeError("PyDESeq2 did not populate DeseqStats.results_df.")
    return results.copy()


def _lfc_shrink_coefficient(dds: Any, contrast: str, group: str) -> str:
    """Resolve the PyDESeq2 0.5 coefficient name for a contrast level."""
    design_matrix = dds.obsm.get("design_matrix")
    columns = list(getattr(design_matrix, "columns", []))
    candidates = [
        f"{contrast}[T.{group}]",
        f"{contrast}_{group}_vs_reference",
    ]
    for candidate in candidates:
        if candidate in columns:
            return candidate

    matching_columns = [
        column
        for column in columns
        if column.startswith(f"{contrast}[T.") and column.endswith("]")
    ]
    if len(matching_columns) == 1:
        return matching_columns[0]

    raise ValueError(
        "Could not identify a PyDESeq2 coefficient for LFC shrinkage. "
        f"Available coefficients: {columns}"
    )


def get_differential_interactions(
    counts: pd.DataFrame,
    design: pd.DataFrame,
    contrast: str,
    group_a: str,
    group_b: str,
    threshold_q: float = 0.05,
    lfc_shrink: bool = False,
) -> pd.DataFrame:
    """Runs DESeq2 on interaction counts."""
    DeseqDataSet, DefaultInference, DeseqStats = _load_pydeseq2()

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
    ds.summary()

    if lfc_shrink:
        ds.lfc_shrink(coeff=_lfc_shrink_coefficient(dds, contrast, group_b))

    df_results = _results_dataframe(ds)

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
        start=lambda df: (
            df.index.str.split(":").str[1].str.split("-").str[0].astype(int)
        ),
        end=lambda df: df.index.str.split(":").str[1].str.split("-").str[1].astype(int),
    )

    return df_results


def differential(
    interaction_files: Sequence[Path | str],
    viewpoint: str,
    design_matrix: Path | str,
    output_prefix: Path | str = "differential_interactions",
    contrast: str = "condition",
    regions_of_interest: Path | str | None = None,
    viewpoint_distance: int | None = None,
    threshold_count: float = 20,
    threshold_q: float = 0.05,
) -> None:
    """Identifies differential interactions between conditions.

    Parses a list of cooler files containg reporter counts from at least two conditions with
    two or more replicates for a single capture probe and outputs differential interaction
    results. Following filtering to ensure that the number of interactions is above the required
    threshold, PyDeseq2 is used to run a compatison after
    fitting a negative binomial model to the interaction counts.The options to filter
    results can be filtered by a minimum mean value (threshold_mean) and/or
    maximum q-value (threshold-q) are also provided.

    Warning:
        Running this on every interaction breaks the model's assumption of
        independence. This is provided as is. For a more statistically sound
        comparison, limit testing to regions of interest.

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
    output_prefix = os.fspath(output_prefix)
    regions_of_interest = (
        os.fspath(regions_of_interest) if regions_of_interest is not None else None
    )

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

    bedgraphs = dict()
    for interaction_file in interaction_files:
        interaction_file = os.fspath(interaction_file)
        file_name = os.path.basename(interaction_file.replace(".hdf5", ""))
        bedgraphs[file_name] = cooler_to_bedgraph(
            clr=f"{interaction_file}::{viewpoint}",
            regions_of_interest=regions_of_interest,
            viewpoint_distance=viewpoint_distance,
        )

    logger.info("Concatenating interactions.")
    # Concatenate bedgraphs
    df_counts = pd.concat(
        [
            bg.assign(
                coord=lambda df: (
                    df["chrom"].astype(str)
                    + ":"
                    + df["start"].astype(str)
                    + "-"
                    + df["end"].astype(str)
                )
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
    if df_counts.empty:
        raise ValueError(
            "No differential interactions found after filtering interactions with "
            f"less than {threshold_count} counts."
        )

    # At the time of writing. PyDeseq2 doese not support multiple comparisons.
    # Therefore, we need to run a separate DESeq2 analysis for each comparison.

    # Get all comparisons
    possible_contrasts = df_design[contrast].unique()
    comparisons = list(itertools.combinations(possible_contrasts, 2))

    # Run comparisons
    for group_a, group_b in comparisons:
        # Filter design matrix
        df_design_sub = df_design.loc[
            lambda df, a=group_a, b=group_b: df[contrast].isin([a, b])
        ]

        # Filter counts
        df_counts_sub = df_counts.loc[:, df_design_sub.index].round().astype(int)

        # Get differential interactions
        logger.info(f"Running comparison: {group_a} vs {group_b}")
        df_results = get_differential_interactions(
            df_counts_sub,
            df_design_sub,
            contrast,
            threshold_q=threshold_q,
            group_a=group_a,
            group_b=group_b,
        )

        # Write result
        df_results.to_csv(
            f"{output_prefix}.{group_a}_vs_{group_b}.csv",
            sep=",",
            index=True,
            header=True,
        )
