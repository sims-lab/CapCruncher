import itertools
import os
import sys
import warnings
from typing import Tuple, Union

import ibis
import numpy as np
import pandas as pd
import pyranges as pr
import pysam
import ray
from loguru import logger
from pybedtools import BedTool, MalformedBedLineError

from capcruncher.api.annotate import BedFileIntersection
from capcruncher.utils import (
    bed_has_name,
    convert_bed_to_dataframe,
    convert_bed_to_pr,
    is_valid_bed,
)

warnings.simplefilter("ignore")


pysam.set_verbosity(0)


def cycle_argument(arg):
    """Allows for the same argument to be stated once but repeated for all files"""

    if len(arg) == 1:
        return itertools.cycle((arg[0],))
    else:
        return arg


def increase_cis_slice_priority(df: pd.DataFrame, score_multiplier: float = 2):
    """
    Prioritizes cis slices by increasing the mapping score.
    """

    df["parent_name"] = df["name"].str.split("|").str[0]

    df_chrom_counts = (
        df[["parent_name", "chrom"]].value_counts().to_frame("slices_per_chrom")
    )
    modal_chrom = (
        df_chrom_counts.groupby("parent_name")["slices_per_chrom"]
        .transform("max")
        .reset_index()
        .set_index("parent_name")["chrom"]
        .to_dict()
    )
    df["fragment_chrom"] = df["parent_name"].map(modal_chrom)
    df["score"] = np.where(
        df["chrom"] == df["fragment_chrom"],
        df["score"] * score_multiplier,
        df["score"] / score_multiplier,
    )

    return df.drop(columns="parent_name")


def remove_duplicates_from_bed(
    bed: pr.PyRanges,
    prioritize_cis_slices: bool = False,
    chroms_to_prioritize: Union[list, np.ndarray] = None,
) -> pr.PyRanges:
    """
    Removes duplicate entries from a PyRanges object.

    Args:
        bed (pr.PyRanges): PyRanges object to be deduplicated.
        prioritize_cis_slices (bool, optional): Prioritize cis slices by increasing the mapping score. Defaults to False.
        chroms_to_prioritize (Union[list, np.ndarray], optional): Chromosomes to prioritize. Defaults to None.

    Returns:
        pr.PyRanges: Deduplicated PyRanges object.
    """

    df = bed.df.rename(columns=lambda col: col.lower()).rename(
        columns={"chromosome": "chrom"}
    )

    # Shuffle the dataframe to randomize the duplicate removal
    df = df.sample(frac=1)

    if prioritize_cis_slices:
        df = increase_cis_slice_priority(df)

    if "score" in df.columns:
        df = df.sort_values(["score"], ascending=False)

    if chroms_to_prioritize:
        df["is_chrom_priority"] = df["chrom"].isin(chroms_to_prioritize).astype(int)
        df = df.sort_values(["score", "is_chrom_priority"], ascending=False).drop(
            columns="is_chrom_priority"
        )

    return (
        df.drop_duplicates(subset="name", keep="first")
        .sort_values(["chrom", "start"])[["chrom", "start", "end", "name"]]
        .rename(columns=lambda col: col.capitalize())
        .rename(columns={"Chrom": "Chromosome"})
        .pipe(pr.PyRanges)
    )


def annotate(
    slices: os.PathLike,
    actions: Tuple = None,
    bed_files: Tuple = None,
    names: Tuple = None,
    overlap_fractions: Tuple = None,
    output: os.PathLike = None,
    duplicates: str = "remove",
    n_cores: int = 1,
    blacklist: str = "",
    prioritize_cis_slices: bool = False,
    priority_chroms: str = "",
    **kwargs,
):
    """
    Annotates a bed file with other bed files using bedtools intersect.

    Whilst bedtools intersect allows for interval names and counts to be used for annotating intervals, this command
    provides the ability to annotate intervals with both interval names and counts at the same time. As the pipeline allows
    for empty bed files, this command has built in support to deal with blank/malformed bed files and will return default N/A values.

    Prior to interval annotation, the bed file to be intersected is validated and duplicate entries/multimapping reads are removed
    to ensure consistent annotations and prevent issues with reporter identification.

    \f
    Args:
     slices (os.PathLike): Input bed file.
     actions (Tuple, optional): Methods to use for annotation. Choose from (get|count). Defaults to None.
     bed_files (Tuple, optional): Bed files to intersect with the bed file to be annotated. Defaults to None.
     names (Tuple, optional): Column names for output tsv file. Defaults to None.
     overlap_fractions (Tuple, optional): Minimum overlap fractions required to call an intersection. Defaults to None.
     output (os.PathLike, optional): Output file path for annotated .tsv file. Defaults to None.
     duplicates (str, optional): Method to deal with multimapping reads/duplicate bed names.
                                 Currently, "remove" is the only supported option. Defaults to "remove".
     n_cores (int, optional): Number of corese to use for intersection. Bed files are intersected in parallel.
                              Defaults to 4.
     invalid_bed_action (str, optional): Action to deal with invalid bed files. Choose from (ignore|error) .These can be ignored by setting to "ignore". Defaults to 'error'.

    Raises:
     NotImplementedError: Only supported option for duplicate bed names is remove.
    """

    with logger.catch():

        logger.info("Validating commandline arguments")
        len_bed_files = len(bed_files)
        if not all([len(arg) == len_bed_files for arg in [actions, names]]):
            raise ValueError(
                "The lengths of the supplied bed files actions and names do not match"
            )

        if slices == "-":
            logger.info("Reading slices from stdin")
            slices = pd.read_csv(sys.stdin, sep="\t", header=None).pipe(pr.PyRanges)

        elif slices.endswith(".bam"):
            logger.info("Converting bam to bed")
            slices = BedTool(slices).bam_to_bed().to_dataframe().pipe(convert_bed_to_pr)

        else:
            slices = pr.PyRanges(slices)

        logger.info("Validating input bed file before annotation")

        if blacklist:
            try:
                logger.info("Removing blacklisted regions from the bed file")
                gr_blacklist = pr.PyRanges(blacklist)
                slices.subtract(gr_blacklist)
            except Exception as e:
                logger.warning(
                    f"Failed to remove blacklisted regions from the bed file. {e}"
                )

        logger.info("Dealing with duplicates in the bed file")

        if duplicates == "remove":
            slices = remove_duplicates_from_bed(
                slices,
                prioritize_cis_slices=prioritize_cis_slices,
                chroms_to_prioritize=priority_chroms.split(",")
                if priority_chroms
                else None,
            )
        else:
            raise NotImplementedError(
                "Only supported option at present is to remove duplicates"
            )

        logger.info("Setting-up intersection(s)")
        ray.init(num_cpus=n_cores, ignore_reinit_error=True, include_dashboard=False)

        # Create a shared object reference to the slices
        slices_ref = ray.put(slices)

        futures = []
        for bed, name, action, fraction in zip(
            bed_files,
            names,
            actions,
            cycle_argument(overlap_fractions),
        ):
            logger.info(f"Setting-up intersection for {bed}")
            bfi = BedFileIntersection.remote(
                bed_a=slices_ref,
                bed_b=bed,
                name=name,
                action=action,
                fraction=fraction,
            )
            futures.append(bfi.intersection.remote())

        # Collate results
        df_annotation = slices.df.set_index("Name")

        while len(futures):
            done_id, futures = ray.wait(futures)

            for ref in done_id:
                ser_new = ray.get(ref)
                df_annotation = df_annotation.join(ser_new, how="left")
                del ser_new

        # Format dataframe for next stage
        df_annotation = (
            df_annotation.reset_index()
            .rename(columns={"Name": "slice_name"})
            .set_index("slice_name")
            .reset_index()
        )

        logger.info("Writing annotations to file.")
        df_annotation.to_parquet(output, compression="snappy")
