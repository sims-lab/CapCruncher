from typing import Union
import pandas as pd
from collections import defaultdict
from itertools import combinations
import xopen
from capcruncher.cli.cli_reporters import cli
import click
import os
from tqdm import tqdm


def count_re_site_combinations(
    groups: pd.core.groupby.GroupBy, column: str = "restriction_fragment", counts: defaultdict = None
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

        for rf1, rf2 in combinations(frag[column], 2):  # Get fragment combinations
            #TODO: Notice a high amount of multicaptures (same oligo) not being removed.
            # Need to track this down but for now will explicitly prevent the same bin appearing twice.
            if not rf1 == rf2:
                rf1, rf2 = sorted([rf1, rf2])  # Sort them to ensure consistency
                counts[rf1, rf2] += 1

    return counts

def count(
    reporters: os.PathLike,
    output: os.PathLike = 'counts.tsv',
    remove_exclusions: bool = False,
    remove_capture: bool = False,
    subsample: int = 0,
):
    """
    Determines the number of captured restriction fragment interactions genome wide.

    Parses a reporter slices tsv and counts the number of unique restriction fragment 
    interaction combinations that occur within each fragment.
    
    Options to ignore unwanted counts e.g. excluded regions or capture fragments are provided. 
    In addition the number of reporter fragments can be subsampled if required.
    
    \f
    Args:
     reporters (os.PathLike): Reporter tsv file path.
     output (os.PathLike, optional): Output file path for interaction counts tsv. Defaults to 'counts.tsv'.
     remove_exclusions (bool, optional): Removes regions marked as capture proximity exclusions before counting. Defaults to False.
     remove_capture (bool, optional): Removes all capture fragments before counting. Defaults to False.
     subsample (int, optional): Subsamples the fragments by the specified fraction. Defaults to 0 i.e. No subsampling.
    """

    try:
        from capcruncher import libcapcruncher
        libcapcruncher.count_restriction_fragment_combinations(reporters, output, 8, int(2e6))
    
    except ImportError as e:

        with xopen.xopen(output, mode="wb", threads=4) as writer:

            # Write output file header.
            header = "\t".join(["bin1_id", "bin2_id", "count"]) + "\n"
            writer.write(header.encode())

            # Open reporter tsv file as a chunked reader (these files can get very large)
            chunksize = 2e6
            df_reporters_iterator = pd.read_csv(reporters, sep="\t", chunksize=chunksize)

            ligated_rf_counts = defaultdict(int)
            for ii, df_reporters in enumerate(df_reporters_iterator):

                print(f'Processing chunk #{ii+1} of {chunksize} slices')

                if remove_exclusions:
                    print('Removing exclusions')
                    # Must only remove exclusions if they are relevant for the capture being examined
                    df_capture = df_reporters.query('capture != "."')

                    # Finds excluded reporters
                    df_reporters_exclusions = df_reporters.query(
                        '(capture == ".") and (exclusion_count > 0)'
                    )

                    # Merge captures with excluded reporters and remove only exclusions
                    # where the excluded region is the same as the capture probe.
                    df_reporters_to_remove = (
                        df_capture.drop_duplicates(["parent_read", "capture"])[
                            ["parent_read", "capture"]
                        ]
                        .merge(
                            df_reporters_exclusions[
                                ["parent_read", "exclusion", "slice_name"]
                            ],
                            on="parent_read",
                        )
                        .query("capture == exclusion")
                    )

                    df_reporters = df_reporters.loc[
                        ~(
                            df_reporters["slice_name"].isin(
                                df_reporters_to_remove["slice_name"]
                            )
                        )
                    ]

                # Remove the capture site
                if remove_capture:
                    print('Removing viewpoints')
                    df_reporters = df_reporters.query('capture != "."')

                # Subsample at the fragment level
                if subsample:
                    print('Subsampling data')
                    if isinstance(subsample, float):
                        subsample_options = {"frac": subsample}
                    elif isinstance(subsample, int):
                        subsample_options = {"n": subsample}

                    # Generate a subsample of fragments and slice these from the reporter dataframe
                    df_reporters = df_reporters[
                        df_reporters["parent_read"].isin(
                            df_reporters["parent_read"].sample(**subsample_options)
                        )
                    ]

                print('Grouping into fragments')
                fragments = df_reporters.groupby("parent_read")

                print('Counting')
                ligated_rf_counts = count_re_site_combinations(
                    fragments, column="restriction_fragment", counts=ligated_rf_counts
                )


            for (rf1, rf2), count in ligated_rf_counts.items():
                line = "\t".join([str(rf1), str(rf2), str(count)]) + "\n"
                writer.write(line.encode())
