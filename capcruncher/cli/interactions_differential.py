"""Identifies differential interactions between conditions."""

import os
import pandas as pd
import itertools
from loguru import logger


import ray
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from capcruncher.api.pileup import cooler_to_bedgraph


@ray.remote
def get_differential_interactions(
    counts: pd.DataFrame,
    design: pd.DataFrame,
    contrast: str,
    threshold_q: float = 0.05,
    lfc_shrink: bool = True,
):
    """Runs DESeq2 on interaction counts."""
    # Create DeseqDataSet
    dds = DeseqDataSet(
        counts=counts.T,
        clinical=design,
        design_factors=contrast,
    )

    # Run DESeq2
    dds.deseq2()

    # Get results
    results = DeseqStats(dds)
    df_results = results.summary()

    if lfc_shrink:
        df_results = results.lfc_shrink()

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
    assert (
        regions_of_interest or viewpoint_distance
    ), "No regions of interest or viewpoint distance provided."

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
            df_counts_sub, df_design_sub, contrast, threshold_q=threshold_q
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


# def differential(
#     union_bedgraph: os.PathLike,
#     capture_name: str,
#     capture_viewpoints: os.PathLike,
#     output_prefix: os.PathLike = "differential",
#     design_matrix: os.PathLike = None,
#     grouping_col: str = "condition",
#     threshold_count: float = 20,
#     threshold_q: float = 0.05,
#     threshold_mean: float = 0,
# ):

#     """
#     Identifies differential interactions between conditions.

#     Parses a union bedgraph containg reporter counts from at least two conditions with
#     two or more replicates for a single capture probe and outputs differential interaction
#     results. Following filtering to ensure that the number of interactions is above the required
#     threshold (--threshold_count), diffxpy is used to run a wald test after
#     fitting a negative binomial model to the interaction counts.The options to filter
#     results can be filtered by a minimum mean value (threshold_mean) and/or
#     maximum q-value (threshold-q) are also provided.

#     Notes:

#      Currently both the capture oligos and the name of the probe being analysed must
#      be provided in order to correctly extract cis interactions.

#      If a N_SAMPLE * METADATA design matrix has not been supplied, the script
#      assumes that the standard replicate naming structure has been followed
#      i.e. SAMPLE_CONDITION_REPLICATE_(1|2).fastq.gz.

#     \f
#     Args:
#         union_bedgraph (os.PathLike): Union bedgraph containg all samples to be compared.
#         capture_name (str): Name of capture probe. MUST match one probe within the supplied oligos.
#         capture_viewpoints (os.PathLike): Capture oligos used for the analysis.
#         output_prefix (os.PathLike, optional): Output prefix for differntial interactions. Defaults to 'differential'.
#         design_matrix (os.PathLike, optional): Design matrix to use for grouping samples. (N_SAMPLES * METADATA). Defaults to None.
#         grouping_col (str, optional): Column to use for grouping. Defaults to 'condition'.
#         threshold_count (float, optional): Minimum number of reported interactions required. Defaults to 20.
#         threshold_q (float, optional): Maximum q-value for output. Defaults to 0.05.
#         threshold_mean (float, optional): Minimum mean value for output. Defaults to 0.
#     """

#     df_bdg = pd.read_csv(union_bedgraph, sep="\t")
#     df_viewpoints = pd.read_csv(
#         capture_viewpoints, sep="\t", names=["chrom", "start", "end", "name"]
#     )

#     #  If design matrix present then use it. Else will assume that the standard format has been followed:
#     #  i.e. NAME_TREATMENT_REPLICATE
#     if design_matrix:
#         df_design = pd.read_csv(design_matrix, sep="\t")
#     else:
#         col_dict = {col: "_".join(col.split("_")[:-1]) for col in df_bdg.columns[3:]}
#         df_design = pd.Series(col_dict).to_frame(grouping_col)

#     # Only cis interactions
#     capture_chrom = get_chromosome_from_name(df_viewpoints, name=capture_name)
#     df_bdg_counts = df_bdg.query(f'chrom == "{capture_chrom}"')

#     # Only counts
#     df_bdg_counts = df_bdg_counts.iloc[:, 3:]

#     # Only with number of interactions > threshold per group in at least 2 replicates
#     df_bdg_counts = (
#         df_bdg_counts.groupby(df_design[grouping_col], axis=1)
#         .apply(lambda df: df[(df >= threshold_count).sum(axis=1) >= 2])
#         .fillna(0.0)
#     )

#     # Run differential testing
#     count_data = df_bdg_counts.transpose().values.astype(np.float64)
#     fragment_names = df_bdg_counts.index.values

#     tests = de.test.pairwise(
#         count_data,
#         grouping=grouping_col,
#         sample_description=df_design,
#         gene_names=fragment_names,
#         test="wald",
#         lazy=False,
#         backend="numpy",
#     )

#     # Go through all of the pairwise tests
#     for g1, g2 in itertools.combinations(tests.groups, 2):
#         df_comparison = tests.summary_pairs(
#             groups0=[
#                 g1,
#             ],
#             groups1=[
#                 g2,
#             ],
#             qval_thres=threshold_q,
#             mean_thres=threshold_mean,
#         )

#         # Extract the fragment coords
#         df_coords = df_bdg.loc[df_comparison["gene"], ["chrom", "start", "end"]]
#         # Merge with test results
#         df_comparison = df_comparison.merge(df_coords, left_on="gene", right_index=True)
#         # Output to tsv
#         (
#             df_comparison.drop(columns="gene")[
#                 ["chrom", "start", "end", "mean", "log2fc", "pval", "qval"]
#             ].to_csv(f"{output_prefix}_{g1}_vs_{g2}.tsv", sep="\t", index=False)
#         )


# # if __name__ == '__main__':
# #     interactions_differential()
