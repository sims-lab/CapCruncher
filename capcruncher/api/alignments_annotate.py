import sys
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Literal

import pandas as pd
import pyranges1 as pr
from loguru import logger

from capcruncher.api.annotate import annotate_intervals, remove_duplicates_from_bed
from capcruncher.utils import (
    convert_bed_to_pr,
    cycle_argument,
    hash_column,
)

warnings.simplefilter("ignore")


def _bam_to_bed_dataframe(bam_path: Path | str) -> pd.DataFrame:
    import pysam

    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_name is None:
                continue

            rows.append(
                {
                    "chrom": read.reference_name,
                    "start": read.reference_start,
                    "end": read.reference_end,
                    "name": read.query_name,
                    "score": read.mapping_quality,
                    "strand": "-" if read.is_reverse else "+",
                }
            )

    return pd.DataFrame.from_records(
        rows, columns=["chrom", "start", "end", "name", "score", "strand"]
    )


def annotate(
    slices: Path | str,
    actions: Sequence[Literal["get", "count"]] | None = None,
    bed_files: Sequence[Path | str] | None = None,
    names: Sequence[str] | None = None,
    overlap_fractions: Sequence[float] | None = None,
    output: Path | str | None = None,
    duplicates: str = "remove",
    n_cores: int = 1,
    blacklist: Path | None = None,
    prioritize_cis_slices: bool = False,
    priority_chroms: list[str] = None,
    **kwargs: Any,
) -> None:
    """
    Annotates a BED-like input with other BED files using PyRanges overlaps.

    Whilst interval overlap tools allow interval names and counts to be used for annotating intervals, this command
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
        actions = tuple(actions or ())
        bed_files = tuple(bed_files or ())
        names = tuple(names or ())
        overlap_fractions = tuple(overlap_fractions or (1e-9,))

        len_bed_files = len(bed_files)
        if not all([len(arg) == len_bed_files for arg in [actions, names]]):
            raise ValueError(
                "The lengths of the supplied bed files actions and names do not match"
            )

        if slices == "-":
            logger.info("Reading slices from stdin")
            slices = convert_bed_to_pr(pd.read_csv(sys.stdin, sep="\t", header=None))

        elif slices.endswith(".bam"):
            logger.info("Converting bam to bed")
            slices = _bam_to_bed_dataframe(slices).pipe(convert_bed_to_pr)

        else:
            slices = convert_bed_to_pr(slices)

        logger.info("Validating input bed file before annotation")

        if blacklist:
            try:
                logger.info("Removing blacklisted regions from the bed file")
                gr_blacklist = pr.read_bed(blacklist)
                slices = slices.subtract_overlaps(gr_blacklist, strand_behavior="ignore")
            except Exception as e:
                logger.warning(
                    f"Failed to remove blacklisted regions from the bed file. {e}"
                )

        logger.info("Dealing with duplicates in the bed file")

        if duplicates == "remove":
            slices: pr.PyRanges = remove_duplicates_from_bed(
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

        for action, bed_file, name, fraction in zip(
            actions, bed_files, names, cycle_argument(overlap_fractions)
        ):
            logger.info(
                f"Performing {name} intersection with {bed_file} using {action} method with {fraction} overlap fraction. {len(slices)} slices to intersect."
            )
            
            slices = annotate_intervals(
                query=slices,
                annotations=bed_file,
                name=name,
                method=action,
                fraction=fraction,
            )
            

        logger.info("Writing annotations to file.")
        df_annotation = slices.rename(columns={"Name": "slice_name"}).assign(
            slice_id=lambda df: hash_column(df.slice_name)
        )
        df_annotation.to_parquet(output, compression="snappy")
