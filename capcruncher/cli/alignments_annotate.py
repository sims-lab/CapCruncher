import itertools
import logging
import sys
import warnings
from typing import Tuple, Union

import numpy as np
from joblib.parallel import Parallel, delayed

warnings.simplefilter("ignore")
import os

import pandas as pd
from capcruncher.tools.annotate import BedIntersection
from capcruncher.utils import bed_has_name, convert_bed_to_dataframe, is_valid_bed
from pybedtools import BedTool, MalformedBedLineError


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
    bed: Union[str, BedTool, pd.DataFrame],
    prioritize_cis_slices: bool = False,
    chroms_to_prioritize: Union[list, np.ndarray] = None,
) -> BedTool:
    """
    Simple removal of duplicated entries from bed file.

    If a "score" field is present a higher scored entry is prioritised.

    Args:
     bed (Union[str, BedTool, pd.DataFrame]): Bed object to deduplicate

    Returns:
     BedTool: BedTool with deduplicated names
    """

    df = convert_bed_to_dataframe(bed).sample(frac=1)

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
        .pipe(BedTool.from_dataframe)
    )


def get_intersection(intersector: BedIntersection):
    return intersector.get_intersection()


def annotate(
    slices: os.PathLike,
    actions: Tuple = None,
    bed_files: Tuple = None,
    names: Tuple = None,
    dtypes: Tuple[bool] = None,
    overlap_fractions: Tuple = None,
    output: os.PathLike = None,
    duplicates: str = "remove",
    n_cores: int = 1,
    invalid_bed_action: str = "error",
    blacklist: str = "",
    prioritize_cis_slices: bool = False,
    priority_chroms: str = "",
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

    logging.info("Validating commandline arguments")
    len_bed_files = len(bed_files)
    if not all([len(arg) == len_bed_files for arg in [actions, names]]):
        raise ValueError(
            "The lengths of the supplied bed files actions and names do not match"
        )

    if slices == "-":
        logging.info("Reading slices from stdin")
        slices = pd.read_csv(sys.stdin, sep="\t", header=None).pipe(
            BedTool.from_dataframe
        )

    elif slices.endswith(".bam"):
        logging.info("Converting bam to bed")
        slices = BedTool(slices).bam_to_bed()

    else:
        slices = BedTool(slices)

    logging.info("Validating input bed file before annotation")
    if not is_valid_bed(slices):
        raise ValueError(f"bed - {slices} is invalid")

    if not bed_has_name(slices):
        raise ValueError(f"bed - {slices} does not have a name column")

    if blacklist:
        logging.info("Removing blacklisted regions from the bed file")
        try:
            slices = slices - BedTool(blacklist)
        except (MalformedBedLineError, FileNotFoundError, IndexError) as e:
            logging.error(f"Blacklist {blacklist} bedfile raised {e}. Ensure it is correctly formatted")


    logging.info("Dealing with duplicates in the bed file")
    # Deal with multimapping reads.
    if duplicates == "remove":
        slices = remove_duplicates_from_bed(slices, 
                                            prioritize_cis_slices=prioritize_cis_slices, 
                                            chroms_to_prioritize=priority_chroms.split(",") if priority_chroms else None)
    else:
        raise NotImplementedError(
            "Only supported option at present is to remove duplicates"
        )

    logging.info("Performing intersection")
    intersections_to_perform = []
    for bed, name, action, fraction, dtype in zip(
        bed_files,
        names,
        actions,
        cycle_argument(overlap_fractions),
        cycle_argument(dtypes),
    ):

        intersections_to_perform.append(
            BedIntersection(
                bed1=slices,
                bed2=bed,
                intersection_name=name,
                intersection_method=action,
                intersection_min_frac=fraction,
                invalid_bed_action=invalid_bed_action,
                dtype=dtype,
            )
        )

    logging.debug(intersections_to_perform)

    intersections_results = Parallel(n_jobs=n_cores)(
        delayed(get_intersection)(intersection)
        for intersection in intersections_to_perform
    )

    logging.info("Merging annotations")

    # Merge intersections with slices
    df_annotation = (
        convert_bed_to_dataframe(slices)
        .rename(columns={"name": "slice_name"})
        .set_index("slice_name")
        .sort_index()
        .join(intersections_results, how="left")
        .reset_index()
        .rename(columns={"index": "slice_name"})
    )

    del intersections_results
    logging.info("Writing annotations to file.")

    df_annotation.loc[
        :, lambda df: df.select_dtypes("number").columns
    ] = df_annotation.select_dtypes("number").astype("float")

    # Export to tsv
    if output.endswith(".tsv"):
        df_annotation.to_csv(output, sep="\t", index=False)
    elif output.endswith(".hdf5"):
        # Need to convert dtypes to ones that are supported
        df_annotation.to_hdf(
            output, key="/annotation", format="table", complib="blosc", complevel=2
        )
    elif output.endswith(".parquet"):
        df_annotation.to_parquet(output, compression="snappy")
