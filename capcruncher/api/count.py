
from typing import Union
import pandas as pd
from collections import defaultdict
import itertools
import xopen
import os
from tqdm import tqdm
import logging


def get_fragment_combinations(df: pd.DataFrame):
    return [sorted(comb) for comb in itertools.combinations(df["restriction_fragment"], 2)]


def subsample_reporters_from_df(df: pd.DataFrame, subsample: float):

    logging.info("Subsampling data")
    if isinstance(subsample, float):
        subsample_options = {"frac": subsample}
    elif isinstance(subsample, int):
        subsample_options = {"n": subsample}

    # Generate a subsample of fragments and slice these from the reporter dataframe
    df_reporters = df[
        df["parent_id"].isin(df["parent_id"].sample(**subsample_options))
    ]

    return df_reporters


def remove_exclusions_from_df(df: pd.DataFrame):
    #TODO: remove this slight dtype hack
    df = df.astype({"viewpoint": str})
    #df.loc[:, "viewpoint"] = df.loc[:, "viewpoint"].astype(str)
    return df.query("viewpoint != exclusion")



def preprocess_reporters_for_counting(df: pd.DataFrame, **kwargs):

    # Need to remove any restiction fragments that are not in the digested genome
    df_reporters = df.query("restriction_fragment != -1")

    if kwargs.get("remove_exclusions"):
        logging.info("Removing excluded regions")
        df_reporters = remove_exclusions_from_df(df_reporters)

    # Remove the capture site
    if kwargs.get("remove_viewpoints"):
        logging.info("Removing viewpoints")
        df_reporters = df_reporters.query('capture_count == 0')

    # Subsample at the fragment level
    if kwargs.get("subsample"):
        df_reporters = subsample_reporters_from_df(df_reporters, kwargs["subsample"])
    
    return df_reporters



def count_re_site_combinations(
    groups: pd.core.groupby.GroupBy,
    column: str = "restriction_fragment",
    counts: defaultdict = None,
) -> defaultdict:
    """
    Counts the number of unique combinations bewteen groups in a column.

    Args:
     df (pd.core.groupby.GroupBy): Aggregated dataframe for processing.
     column (str, optional): Column to examine for unique combinations per group. Defaults to "restriction_fragment".
     counts (defaultdict, optional): defaultdict(int) containing previous counts. Defaults to None.

    Returns:
     defaultdict: defaultdict(int) containing the count of unique interactions.
    """

    if not counts:
        counts = defaultdict(int)  # Store counts in a default dict

    # For each set of ligated fragments
    for ii, (group_name, frag) in enumerate(tqdm(groups)):

        for rf1, rf2 in itertools.combinations(frag[column], 2):  # Get fragment combinations
            # TODO: Notice a high amount of multicaptures (same oligo) not being removed.
            # Need to track this down but for now will explicitly prevent the same bin appearing twice.
            if not rf1 == rf2:
                rf1, rf2 = sorted([rf1, rf2])  # Sort them to ensure consistency
                counts[rf1, rf2] += 1

    return counts


def get_counts_from_tsv_by_batch(reporters: os.PathLike, chunksize: int, **kwargs):

    df_reporters_iterator = pd.read_csv(reporters, sep="\t", chunksize=chunksize)

    ligated_rf_counts = defaultdict(int)
    for ii, df_reporters in enumerate(df_reporters_iterator):

        logging.info(f"Processing chunk #{ii+1} of {chunksize} slices")

        reporters = preprocess_reporters_for_counting(df_reporters, **kwargs)
        fragments = reporters.groupby("parent_id")

        logging.info("Counting")
        ligated_rf_counts = count_re_site_combinations(
            fragments, column="restriction_fragment", counts=ligated_rf_counts
        )

    return ligated_rf_counts


def get_counts_from_tsv(reporters: os.PathLike, **kwargs):

    df_reporters = pd.read_csv(reporters, sep="\t")

    reporters = preprocess_reporters_for_counting(df_reporters, **kwargs)
    fragments = reporters.groupby("parent_id")

    logging.info("Counting")
    ligated_rf_counts = count_re_site_combinations(
        fragments, column="restriction_fragment"
    )

    return ligated_rf_counts