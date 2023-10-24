import os
import sys
import warnings
from typing import Tuple

import pandas as pd
import pyranges as pr
from loguru import logger
from pybedtools import BedTool

from capcruncher.api.annotate import remove_duplicates_from_bed, BedIntersector
from capcruncher.utils import (
    convert_bed_to_pr,
    cycle_argument
)

warnings.simplefilter("ignore")

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
        
        for action, bed_file, name, fraction in zip(actions, bed_files, names, cycle_argument(overlap_fractions)):
            logger.info(f"Performing {name} intersection with {bed_file} using {action} method with {fraction} overlap fraction. {len(slices)} slices to intersect.")
            slices = BedIntersector(
                bed_a=slices,
                bed_b=bed_file,
                name=name,
                fraction=fraction,
            ).get_intersection(method=action)
            
        

        logger.info("Writing annotations to file.")
        df_annotation = slices.df.rename(columns={"Name": "slice_name"}).assign(slice_id=lambda df: df.slice_name.astype("category"))
        df_annotation.to_parquet(output, compression="snappy")
